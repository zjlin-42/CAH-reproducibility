function zoomin3D(tau, N, direction, p, q)
% ZOOMIN3D - Study the hysteresis under rotation of rough surface, in a zoom-in domain
%
%   zoomin3D(tau, N, direction, p, q)
%
% Inputs:
%   tau: time step, positive
%   N: size of meshed grids, positive even number
%   direction: 'left' or 'right'
%   p, q: integers, atan(p/q) as the angle of rotation


% Constants and Parameters
maxit = 40000;

R = 1.5;

epsilon = 10^-1; % roughness scalar
l_x = 4 * pi * epsilon; 
l_z = 6 * epsilon;
r = 3 * epsilon;

tol=10^-7*N^3; % tolerance
% sigma_LV = 1; % liquid–vapor surface energy density
theta = pi/3; % Young's angle

time_step = tau;
tau2 = (pi*cos(theta)/(pi-2*theta))^2*tau;

gcd_pq = gcd(p,q);
p = p/gcd_pq;
q = q/gcd_pq;
rotation = atan2(p,q);
l_y = 2*pi*epsilon*sqrt(p^2+q^2)/2;

direction = lower(string(direction));
switch direction
    case "left"
        k=3; %inverse of slope of [0 R] to [upper_x upper_z], i.e. k=(0-x)/(R-z)
    case "right"
        k=-3; %inverse of slope of [0 R] to [upper_x upper_z], i.e. k=(0-x)/(R-z)
    otherwise
        error('Invalid direction input');
end


% centroid point (center_x, center_y=0, center_z=-R) and mid point (on \partial C_far)
mid_z=-R+2/3*r; %height of C_near\C_far

center_x=-(2*R-r)*k; % initialzation of center_x; used to shift the cube

% Generate grid points
x = linspace(-2 * l_x, 2 * l_x, 2 * N + 1);
y = linspace(-l_y, l_y, N + 1);
z = linspace(-l_z, l_z, N + 1);

[X, Y, Z] = meshgrid(x(1:end-1), y(1:end-1), z(1:end-1));
dx = x(2) - x(1);
dy = y(2) - y(1);
dz = z(2) - z(1);


% Set heat kernel
f = 1/(4*pi*tau)^(3/2)*exp(-(X.^2+Y.^2+Z.^2)/4/tau); 
g = 1/(4*pi*tau2)^(3/2)*exp(-(X.^2+Y.^2+Z.^2)/4/tau2);
Fhat = fftn(f);
Ghat = fftn(g);

% initial area of liquid D1, initial area of vapor D2, area of solid D3
D_left = double(X<=0);
D_right = ~D_left;

% set rotated solid domain
phi = @(X) -r+epsilon*sin(X{1}/epsilon).*sin(X{2}/epsilon);
XY_rotated = @(X) {X{1} * cos(rotation) - X{2} * sin(rotation), X{1} * sin(rotation) + X{2} * cos(rotation)};
D3 = (double(Z <= phi(XY_rotated({X+center_x+l_x, Y}))) & D_left) + (double(Z <= phi(XY_rotated({l_x-X+center_x, Y}))) & D_right);

D4 = double(Z>=r); % Imaginary upper solid
D5 = (Z>=2/3*r); % Far region

D1 = double(X<-l_x | X>l_x) & ~D3 & ~D4; % Liquid

D2 = ~D1 & ~D3 & ~D4; % Vapor

% Initialize data to be stored
T=0; % total time
D_star=D1; % for modified algorithm convergence
K=zeros(maxit,1); % record all slopes
record_x = zeros(maxit,1); % record of middle_x
record_x_mid = zeros(maxit,1); % record of top_x
convergence = false;

for i=1:maxit
    interface = circshift(D2,[0 -1 0]).*D1; % only consider regular interface
    [~,preshift_x,~]=ind2sub(size(interface),find(interface));
    preshift_x=sum(preshift_x,'all')/sum(interface,'all');
    preshift_x=floor(preshift_x-N/2-1);
    
    if isempty(preshift_x) || sum(interface,'all')==0 % In consideration of the case of no interface
        % Assign a default value (in this case, assigning 0)
        preshift_x = 0;
    end
    
    % update all data after shifting

    center_x = center_x + preshift_x*dx; % is actually the average x of the interface
    record_x(i) = center_x;

    D3 = (double(Z <= phi(XY_rotated({X+center_x+l_x, Y}))) & D_left) + (double(Z <= phi(XY_rotated({l_x-X+center_x, Y}))) & D_right);

    D1 = (circshift(D1 & D_left & double(X>=-2*l_x+dx*preshift_x-0.5*dx),[0 -preshift_x 0])) + (double(X<-2*l_x-dx*preshift_x-0.5*dx) & ~D3 & ~D4 & D_left) + (circshift(D1 & D_right & double(X<2*l_x-dx*preshift_x-0.5*dx),[0 preshift_x 0])) + (double(X>=2*l_x+dx*preshift_x-0.5*dx) & ~D3 & ~D4 & D_right);
    D2 = ~D1 & ~D3 & ~D4;

    D_star = (circshift(D_star & D_left & double(X>=-2*l_x+dx*preshift_x-0.5*dx),[0 -preshift_x 0])) + (double(X<-2*l_x-dx*preshift_x-0.5*dx) & ~D3 & ~D4 & D_left) + (circshift(D_star & D_right & double(X<2*l_x-dx*preshift_x-0.5*dx),[0 preshift_x 0])) + (double(X>=2*l_x+dx*preshift_x-0.5*dx) & ~D3 & ~D4 & D_right);

    % construct the upper imaginary domain and record value of k
    mid_interface = circshift(D2,[0 -1 0]) & D1 & ~D5; % have to update the interface after shifting
    mid_interface = circshift(mid_interface,[0 0 1]) & D5;
    [mid_y,mid_x,~]=ind2sub(size(mid_interface),find(mid_interface)); % note that the data is local location
    [~, sort_order] = sort(mid_y);
    mid_x = mid_x(sort_order);
    mid_x = center_x + (mid_x-N/2-1)*dx; % now mid_x is stored as global location
    mean_mid_x = mean(mid_x);
    record_x_mid(i) = mean_mid_x;

    k=(0-mean_mid_x)/(R-mid_z); %inverse of the slope on average
    K(i)=k;

    if length(mid_x) == N
        % Compute coefficient for truncated U=a0*z+\sum_{n=1}^{N} [an*cos(2*pi*n*y/L)+bn*sin(2*pi*n*y/L)]*sinh(2*pi*n*z/L)
        [a0,an,bn]=solve_Laplace(N,2*l_y/sqrt(1+k^2),2*R-2/3*r,mid_x-mean_mid_x);
        temp_imagnary = imaginary_interface(Y,Z-R,R,l_y,k,a0,an,bn);
        % update portion of D1,D2 over D5\D4
        D51 = (D5 & (1-D4) & D_left & double(X+l_x+center_x < temp_imagnary)) + (D5 & D_right & double(-X+center_x+l_x < temp_imagnary));
        D52 = ~D51 & D5 & (1-D4);
    else
        error("Irregular mid interface detected.")
    end
    
    T=T+tau;

    % compute convolution and take thresholding
    temp=D1;
    % update parameters based on theta_i
    theta_i = pi/2-atan(k);
    tau3 = (pi*cos(theta_i)/(pi-2*theta_i))^2*tau;
    h = 1/(4*pi*tau3)^(3/2)*exp(-(X.^2+Y.^2+Z.^2)/4/tau3);
    % convolution
    % factor dx*dy*dz omitted since the threshold is zero
    u = 1/sqrt(tau)*ifftshift(ifftn(Fhat.*fftn((D2-D1).*(1-D5)+D52-D51))) - cos(theta)/sqrt(tau2)*ifftshift(ifftn(Ghat.*fftn(D3))) - cos(theta_i)/sqrt(tau3)*ifftshift(ifftn(fftn(h).*fftn(D4)));
    
    % delta is fixed 0 without volume constrant; delta should change with each itration under volume constrant
    delta=0; 
    D1 = double(u<delta) & ~D3 & ~D4;
    D2 = ~D1 & ~D3 & ~D4;


    % stopping criteria 
    if sum(abs(D1-temp),"all")<=tol % check difference with last diffusion step
        if sum(abs(D1-D_star),"all")>=tol % check difference with last time step change
            tau=tau/2;
            tau2=tau2/2;
            f = 1/(4*pi*tau)^(3/2)*exp(-(X.^2+Y.^2+Z.^2)/4/tau);
            g = 1/(4*pi*tau2)^(3/2)*exp(-(X.^2+Y.^2+Z.^2)/4/tau2); 
            Fhat = fftn(f);
            Ghat = fftn(g);
            D_star=D1;
        else
            disp("iteration step:" + string(i))
            disp("total time of evolution: " + string(T))
            convergence = true;
            break
        end
    end
end

% store data and figure

fig = figure;

isoD1=isosurface(X(:,1:N,1:N),Y(:,1:N,1:N),Z(:,1:N,1:N),D1(:,1:N,1:N));
p1 = patch(isoD1);
set(p1,'EdgeColor','none');
view(3);
set(p1,'FaceColor',[0 0.5 0.9]);
set(p1,"FaceAlpha",1);

hold on

camlight;
lighting gouraud;

isoD51=isosurface(X(:,1:N,1:N),Y(:,1:N,1:N),Z(:,1:N,1:N),D51(:,1:N,1:N));
p51 = patch(isoD51);
set(p51,'EdgeColor','none');
set(p51,'FaceColor',[0 0.5 0.9]);
set(p51,"FaceAlpha",0.7);

isoD3=isosurface(X(:,1:N,1:N),Y(:,1:N,1:N),Z(:,1:N,1:N),D3(:,1:N,1:N));
p3 = patch(isoD3);
set(p3,'EdgeColor','none');
set(p3,'FaceColor',[1 1 1]);
set(p3,"FaceAlpha",0.3);


% Suggested graphical setting

% axis normal;
% view(60,30);
% axis off;
% h = findobj(gca,'Type','light');
% h.Color = [1 1 1];
% h.Style = 'infinite';
% h.Position = [1 -1 1];
% material shiny;


final_contact_angle = pi/2+atan(k);


foothill = D1.*circshift(D3,[0 0 1]);%.*double(-r-eps<=Z).*double(Z<=-r+eps);
foothill = sum(foothill,3);
foothill = foothill(:,1:N);
fig_foothill = figure;
hold on
contour(X(:,1:N,1),Y(:,1:N,1),foothill,1,'LineWidth',2,'LineColor','red')
sum_D3 = sum(D3(:,1:N,:),3);
contour(X(:,1:N,1),Y(:,1:N,1),sum_D3(:,1:N),'LineWidth',1,'LineColor','blue')

%threshold of u
u=(u>0);

% Save fig directly may cause sleeping for serval directions
% save('3D_'+string(direction)+'_'+string(p)+'_'+string(q)+'_'+regexprep(sprintf('%.0e',time_step),'e([-+])0*(\d+)','e$1$2')+'_zoomin.mat','K','mean_mid_x','center_x','final_contact_angle','record_x_mid','record_x','D1','D51','D3','u','fig','fig_foothill','foothill','-v7.3')

% Save iso's instead
% save('3D_'+string(direction)+'_'+string(p)+'_'+string(q)+'_'+regexprep(sprintf('%.0e',time_step),'e([-+])0*(\d+)','e$1$2')+'_zoomin.mat','K','mean_mid_x','center_x','final_contact_angle','record_x_mid','record_x','D1','D51','D3','u','isoD1','isoD51','isoD3','fig_foothill','foothill','-v7.3')

% Save only the necessary
save('3D_'+string(direction)+'_'+string(p)+'_'+string(q)+'_'+regexprep(sprintf('%.0e',time_step),'e([-+])0*(\d+)','e$1$2')+'_zoomin.mat','K','mean_mid_x','center_x','final_contact_angle','record_x_mid','record_x','D1','D51','D3','u','fig_foothill','foothill','convergence','-v7.3')


close all;

% warning for convergence failure
if ~convergence
    disp("Warning: task might not converge.")
end

end








function [a0,an,bn]=solve_Laplace(N,L,H,f)
% Solve Laplace equation in rectangle domain ([0,L]x[0,H]) and find uy (u partial y) at y=H
% Input:
% N: number of grids
% L: period of x axis (peridic boundary condition)
% H: range y from 0 to H
% f: boundary condition at y = H 
% Other B.C.: u(x,0)=0 and u(0,y)=u(L,y)
% Output: Nx1 vector uy

f = reshape(f, 1, []); % turn f into a row vector
F = fft(f); %DFT of f

% Compute coefficient truncated u:=a0*y+\sum_{n=1}^{N} [an*cos(2*pi*n*x/L)+bn*sin(2*pi*n*x/L)]*sinh(2pi*n*y/L)
a0 = 1/N/H*F(1);
an = 2./sinh(2*pi*(1:N-1)*H/L)/N.*real(F(2:N));
bn = -2./sinh(2*pi*(1:N-1)*H/L)/N.*imag(F(2:N));


end


function [X]=imaginary_interface(Y,Z,R,l_y,k,a0,an,bn)
% Reconstruct epsilon*u(y,z)=U(1/sqrt(1+k^2)*(y+l_y),R-z)
% = a0*(R-z)+\sum_{n=1}^{N} [an*cos(2*pi*n*1/sqrt(1+k^2)*(y+l_y)/L)+bn*sin(2*pi*n*1/sqrt(1+k^2)*(y+l_y)/L)]*sinh(2*pi*n*(R-z)/L)
% and compute x(y,z)=k(z-R)+epsilon*u(y,x).
if ~isequal(size(Y),size(Z))
    error("Arrays have incompatible sizes for this operation.")
end

X=zeros(size(Y));
X=X+k*(Z-R)+a0*(R-Z);
for i=1:find(abs(an)>eps,1,'last')
    X=X+an(i).*cos(pi*i*(Y+l_y)/l_y).*sinh(pi*i*(R-Z)/l_y*sqrt(1+k^2));
end
for i=1:find(abs(bn)>eps,1,'last')
    X=X+bn(i).*sin(pi*i*(Y+l_y)/l_y).*sinh(pi*i*(R-Z)/l_y*sqrt(1+k^2));
end

end
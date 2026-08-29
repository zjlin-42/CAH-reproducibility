%% Demo for directional contact angle hysteresis

epsilon = 1e-2;
N = 64;
p = 6;
q = 6;

% Reduce the direction to relatively prime integers
gcd_pq = gcd(p,q);
p = p/gcd_pq;
q = q/gcd_pq;

% Receding case
zoomin3D(epsilon, N, "right", p, q);
filename_rec = sprintf('3D_right_%d_%d_1e-3_zoomin.mat', p, q);
data_rec = load(filename_rec);
theta_rec = data_rec.final_contact_angle;

% Advancing case
zoomin3D(epsilon, N, "left", p, q);
filename_adv = sprintf('3D_left_%d_%d_1e-3_zoomin.mat', p, q);
data_adv = load(filename_adv);
theta_adv = data_adv.final_contact_angle;

% Display the results
fprintf('Direction: atan2(%d,%d) = %.6f rad\n', p, q, atan2(p,q));
fprintf('Receding contact angle:  %.6f rad\n', theta_rec);
fprintf('Advancing contact angle: %.6f rad\n', theta_adv);
fprintf('CAH interval: [%.6f, %.6f] rad\n', theta_rec, theta_adv);
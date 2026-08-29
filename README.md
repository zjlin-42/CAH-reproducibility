# Reproducibility Files for the Contact Angle Hysteresis Simulations

This repository contains MATLAB code for a reduced-resolution demonstration of the directional contact angle hysteresis (CAH) calculation.
The demonstration is based on the periodic rough surface considered in the manuscript.

## Files

- `zoomin3D.m`: Main simulation function for computing a receding or advancing contact angle at a selected contact-line orientation.
- `demo.m`: Runs one representative directional CAH calculation.

## Requirements

MATLAB is required. No external input data files are needed.

## Usage

Download or clone the repository, open the repository folder in MATLAB, and set it as the current folder. Then run

```matlab
demo
```

The demonstration uses

- initial time step `tau = 1e-2`;
- spatial resolution `N = 64`; and
- contact-line orientation
  `arg(k) = atan2(1,1) = pi/4`.

The script runs two simulations:

- `"right"` for the receding contact angle;
- `"left"` for the advancing contact angle.

After both simulations finish, the script displays the receding angle, advancing angle, CAH interval, and interval width in radians.
The following output files are saved in the current MATLAB folder:

```text
3D_right_1_1_1e-2_zoomin.mat
3D_left_1_1_1e-2_zoomin.mat
```

This reduced-resolution example is intended to demonstrate the computational workflow. 
The numerical experiments reported in the manuscript were performed at a higher spatial resolution.

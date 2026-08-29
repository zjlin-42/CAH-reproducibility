# Reproducibility Files for the Contact Angle Hysteresis Simulations

This repository contains the MATLAB code used to compute the directional
contact angle hysteresis (CAH) for the periodic rough surface considered in
the manuscript.

## Files

- `zoomin3D.m`: Main simulation code for computing the contact angle.
- `demo.m`: A quick demonstration of the directional CAH calculation.

## Usage

Run `demo.m` in MATLAB. The demo runs a representative case with time step
`tau = 1e-2` and spatial resolution `N = 64`, and displays the receding and
advancing contact angles for a selected direction.

The simulations used in the manuscript were performed with higher spatial
resolution.

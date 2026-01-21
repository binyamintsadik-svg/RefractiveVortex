# Refractive Vortex: Resolving the Hubble Tension

## Overview
This repository contains the codebase and data analysis pipeline for the paper **"Resolving the Hubble Tension via Refractive Lensing: A Cubic-Root Dispersion Law."**

We demonstrate that the Hubble Tension can be resolved by modeling the local Large Scale Structure (LSS) as a refractive medium with a frequency-dependent refractive index scaling as $\alpha = 1/3$. This model reconciles the local Supernova measurements (Pantheon+, SH0ES) with the CMB acoustic scale (Planck) without modifying the standard background expansion ($H_0 \approx 67.5$ km/s/Mpc).

## Repository Structure

```text
/
├── physics_model.py        # Core refractive engine (NFW potential, integration)
├── main_analysis.py        # Grid scans, likelihood minimization, and Monte Carlo
├── figure1_plot.py         # Generates the "Holy Trinity" Hubble Diagram
└── requirements.txt        # Python dependencies
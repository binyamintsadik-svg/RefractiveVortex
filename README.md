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


## Reproducing paper figures

### Figure 1 (hybrid decomposition)

This repository includes `figure1_inputs.npz` which contains the precomputed arrays needed to reproduce Figure 1 exactly.

Run:

python plot_hybrid_decomposition.py --inputs figure1_inputs.npz

Outputs are written to `figures/`.

### Appendix Figure A1 (teardrop geometry)

Run:

python plot_teardrop_geometry.py --rs1 830 --rs2 2000

Outputs are written to `figures/`.

### Production MCMC statistics

See `production_run_stats.md` for the posterior summaries from the 12k production run.

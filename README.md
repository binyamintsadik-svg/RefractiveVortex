# Refractive Vortex: Resolving the Hubble Tension

## Overview
This repository contains the codebase and data analysis pipeline for the paper **"Resolving the Hubble Tension via Refractive Lensing: A Cubic-Root Dispersion Law."**

We demonstrate that the Hubble Tension can be resolved by modeling the local Large Scale Structure (LSS) as a refractive medium with a frequency-dependent refractive index scaling as $\alpha = 1/3$. This model reconciles the local Supernova measurements (Pantheon+, SH0ES) with the CMB acoustic scale (Planck) without modifying the standard background expansion ($H_0 \approx 67.5$ km/s/Mpc).

## Repository Structure

```text
/
├── README.md                <-- Update this with the new Abstract & Headline Numbers
├── paper/
│   └── cosmological_lensing.pdf  <-- The final compiled PDF
├── code/
│   ├── plotting/
│   │   ├── plot_hybrid_decomposition.py  <-- NEW (Fig 1)
│   │   └── plot_teardrop_geometry.py     <-- NEW (Appendix A)
│   ├── analysis/
│   │   └── mcmc_production_run.py        <-- The 12k Step Runner
│   └── utils/
│       └── physics_model.py              <-- The core integral logic
├── data/
│   ├── figure1_inputs.npz                <-- The Pantheon+ data subset
│   └── production_run_summary.txt        <-- NEW (The stats)
└── notebooks/
    └── jwst_age_paradox.ipynb            <-- Optional: The z=14 vs z=10 math

"""plot_hybrid_decomposition.py

Generate the paper's Figure 1: hybrid decomposition stacked panels (shift + residuals).

This script is intentionally self-contained and "paper reproducibility" grade.

Inputs
- figure1_inputs.npz (committed): precomputed arrays exported from the analysis pipeline.

Expected NPZ keys (minimum)
- z_sorted: sorted redshifts for SNe
- res_sn: SN residuals (mag) relative to baseline model used in the pipeline
- err_sorted: SN residual uncertainties (mag)
- z_line: redshift grid for smooth curves
- res_shoes_line: reference line residual (e.g., SH0ES)
- res_ref_line: model line residual (e.g., CLE / refractive)

Optional NPZ keys (if present)
- bao_z, bao_res: BAO points for overlay

Outputs
- figures/figure1_hybrid_decomposition.png
- figures/figure1_hybrid_decomposition.pdf

Usage
python plot_hybrid_decomposition.py --inputs figure1_inputs.npz

Notes
- This script only plots. It does not refit parameters.
- If your pipeline changes variable naming, update the NPZ exporter to match the keys above.
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def _load_inputs(npz_path):
    dat = np.load(npz_path)
    keys = list(dat.keys())

    required = ['z_sorted', 'res_sn', 'err_sorted', 'z_line', 'res_shoes_line', 'res_ref_line']
    missing = [k for k in required if k not in keys]
    if len(missing) > 0:
        raise KeyError('Missing required NPZ keys: ' + ', '.join(missing) + '. Found keys: ' + ', '.join(keys))

    out = {k: dat[k] for k in required}

    # Optional overlays
    if 'bao_z' in keys and 'bao_res' in keys:
        out['bao_z'] = dat['bao_z']
        out['bao_res'] = dat['bao_res']

    return out


def make_figure(inputs, out_png, out_pdf):
    z_sorted = inputs['z_sorted']
    res_sn = inputs['res_sn']
    err_sorted = inputs['err_sorted']

    z_line = inputs['z_line']
    res_shoes_line = inputs['res_shoes_line']
    res_ref_line = inputs['res_ref_line']

    # Two stacked panels
    fig = plt.figure(figsize=(8.0, 7.5))

    ax1 = plt.subplot(2, 1, 1)
    ax2 = plt.subplot(2, 1, 2, sharex=ax1)

    # Top: shift / comparison lines
    ax1.plot(z_line, res_shoes_line, color='tab:orange', lw=2.0, label='Reference (SH0ES)')
    ax1.plot(z_line, res_ref_line, color='tab:blue', lw=2.0, label='Model (CLE / refractive)')

    # Optional BAO overlay
    if 'bao_z' in inputs:
        ax1.scatter(inputs['bao_z'], inputs['bao_res'], s=24, color='k', alpha=0.8, label='BAO (overlay)')

    ax1.set_ylabel('Residual (mag)')
    ax1.set_title('Hybrid Decomposition (Shift + Residuals)')
    ax1.grid(True, alpha=0.25)
    ax1.legend(frameon=False, loc='best')

    # Bottom: SN residual scatter
    ax2.errorbar(z_sorted, res_sn, yerr=err_sorted, fmt='.', ms=3.5, color='k', alpha=0.65, ecolor='0.7', elinewidth=0.6)
    ax2.axhline(0.0, color='0.2', lw=1.0)
    ax2.set_xlabel('Redshift z')
    ax2.set_ylabel('SN residuals (mag)')
    ax2.grid(True, alpha=0.25)

    plt.tight_layout()

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220, bbox_inches='tight')
    fig.savefig(out_pdf, bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--inputs', type=str, default='figure1_inputs.npz')
    ap.add_argument('--outdir', type=str, default='figures')
    args = ap.parse_args()

    npz_path = Path(args.inputs)
    outdir = Path(args.outdir)

    inputs = _load_inputs(npz_path)

    out_png = outdir / 'figure1_hybrid_decomposition.png'
    out_pdf = outdir / 'figure1_hybrid_decomposition.pdf'

    make_figure(inputs, out_png, out_pdf)


if __name__ == '__main__':
    main()

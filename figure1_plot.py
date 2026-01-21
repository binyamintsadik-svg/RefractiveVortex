"""figure1_plot.py

Creates the publication Figure 1 (Hubble residual diagram) and saves PNG/PDF.

This script expects precomputed arrays saved to NPZ (for exact reproduction from the notebook).
You can export those arrays from your analysis pipeline.

Expected NPZ keys:
- z_sorted, res_sn, err_sorted
- bao_z, bao_res
- z_line, res_shoes_line, res_ref_line

Outputs:
- figure1_grand_unified_polished.png
- figure1_grand_unified_polished.pdf
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def main(npz_in='figure1_inputs.npz'):
    d = np.load(npz_in)
    z_sorted = d['z_sorted']
    res_sn = d['res_sn']
    err_sorted = d['err_sorted']
    bao_z = d['bao_z']
    bao_res = d['bao_res']
    z_line = d['z_line']
    res_shoes_line = d['res_shoes_line']
    res_ref_line = d['res_ref_line']

    plt.rcParams.update({
        'font.size': 13,
        'axes.titlesize': 16,
        'axes.labelsize': 15,
        'legend.fontsize': 11,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'figure.dpi': 120,
        'savefig.dpi': 450,
        'font.family': 'serif',
        'mathtext.fontset': 'cm'
    })

    fig, ax_main = plt.subplots(figsize=(12.5, 6.8))

    ax_main.errorbar(z_sorted, res_sn, yerr=err_sorted, fmt='.', color='0.35', alpha=0.18,
                     markersize=2.2, linewidth=0.45, label='Pantheon+SH0ES SNe Ia')
    ax_main.scatter(bao_z, bao_res, marker='D', s=70, color='dodgerblue', edgecolor='black', linewidth=0.6,
                    label='BAO')
    ax_main.scatter([1100.0], [0.0], marker='*', s=320, color='gold', edgecolor='black', linewidth=0.8,
                    label='CMB anchor (z=1100)')

    ax_main.plot(z_line, np.zeros_like(z_line), linestyle='--', color='black', linewidth=1.8,
                label='Planck LCDM baseline (H0=67.4)')
    ax_main.plot(z_line, res_shoes_line, linestyle=':', color='red', linewidth=2.2,
                label='SH0ES-like local slope (H0=73.0)')
    ax_main.plot(z_line, res_ref_line, color='green', linewidth=3.2,
                label='Refractive vortex model')

    ax_main.set_xscale('log')
    ax_main.set_xlim(0.01, 1100.0)
    ax_main.set_xlabel('Redshift z')
    ax_main.set_ylabel('Distance modulus residual  mu_obs - mu_Planck  mag')
    ax_main.set_title('Refractive Reconciliation: SNe Ia + BAO + CMB on a unified residual Hubble diagram')
    ax_main.axhline(0.0, color='black', linewidth=0.8, alpha=0.5)
    ax_main.legend(loc='upper right', frameon=True, framealpha=0.92)

    ax_in = inset_axes(ax_main, width='42%', height='48%', loc='lower left', borderpad=2.0)
    mask_low = z_sorted < 0.15
    ax_in.errorbar(z_sorted[mask_low], res_sn[mask_low], yerr=err_sorted[mask_low], fmt='.',
                   color='0.35', alpha=0.25, markersize=3.0, linewidth=0.6)

    mask_line_low = z_line <= 0.15
    ax_in.plot(z_line[mask_line_low], np.zeros_like(z_line[mask_line_low]), linestyle='--', color='black', linewidth=1.2)
    ax_in.plot(z_line[mask_line_low], res_shoes_line[mask_line_low], linestyle=':', color='red', linewidth=1.8)
    ax_in.plot(z_line[mask_line_low], res_ref_line[mask_line_low], color='green', linewidth=2.6)

    ax_in.set_xlim(0.01, 0.15)
    low_vals = res_sn[mask_low]
    ax_in.set_ylim(float(np.nanpercentile(low_vals, 1)) - 0.05, float(np.nanpercentile(low_vals, 99)) + 0.05)
    ax_in.set_title('Zoom: z<0.15', fontsize=12)
    ax_in.grid(True, alpha=0.2)

    fig.subplots_adjust(left=0.08, right=0.98, top=0.90, bottom=0.12)

    fig.savefig('figure1_grand_unified_polished.png', bbox_inches='tight')
    fig.savefig('figure1_grand_unified_polished.pdf', bbox_inches='tight')


if __name__ == '__main__':
    main()

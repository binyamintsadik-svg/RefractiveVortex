"""plot_teardrop_geometry.py

Generate Appendix Figure A1: compare the phi(z) profile ("teardrop" geometry) for two scale radii.

This script is meant to visualize how the normalized profile varies with r_s.

Assumptions
- Uses the same profile form as physics_model.py: phi(x) = -ln(1+x)/x, with x = r/r_s.
- Uses flat LCDM to map redshift z -> comoving distance r(z) for a qualitative geometry plot.

Outputs
- figures/figureA1_teardrop_geometry.png
- figures/figureA1_teardrop_geometry.pdf

Usage
python plot_teardrop_geometry.py --rs1 830 --rs2 2000

Note
- This is a visualization tool, not a parameter-fitting script.
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def Ez(z, Om=0.3, Ol=0.7):
    return np.sqrt(Om * (1.0 + z) ** 3 + Ol)


def comoving_distance_Mpc(z_grid, H0=67.4, Om=0.3):
    # Simple numerical integral for D_C(z) = c/H0 * 
    c_km_s = 299792.458
    Ol = 1.0 - Om
    integrand = 1.0 / Ez(z_grid, Om=Om, Ol=Ol)
    dz = z_grid[1] - z_grid[0]
    cum = np.cumsum((integrand[:-1] + integrand[1:]) * 0.5 * dz)
    cum = np.insert(cum, 0, 0.0)
    return (c_km_s / H0) * cum


def phi_profile(x):
    # Stable evaluation near x~0
    x = np.asarray(x)
    out = np.zeros_like(x, dtype=float)
    small = x < 1e-6
    out[small] = -1.0
    xs = x[~small]
    out[~small] = -np.log(1.0 + xs) / xs
    return out


def make_figure(rs1, rs2, out_png, out_pdf):
    z_grid = np.linspace(0.0, 3.0, 2500)
    r_grid = comoving_distance_Mpc(z_grid)

    x1 = r_grid / float(rs1)
    x2 = r_grid / float(rs2)

    phi1 = phi_profile(x1)
    phi2 = phi_profile(x2)

    fig = plt.figure(figsize=(8.0, 4.8))
    plt.plot(z_grid, phi1, lw=2.2, label='r_s = ' + str(rs1) + ' Mpc')
    plt.plot(z_grid, phi2, lw=2.2, label='r_s = ' + str(rs2) + ' Mpc')

    # Mark the approximate z for D_C=830 Mpc (visual cue only)
    # We find z where r(z) is closest to 830 Mpc.
    idx = int(np.argmin(np.abs(r_grid - 830.0)))
    z_830 = float(z_grid[idx])
    plt.axvline(z_830, color='0.2', lw=1.2, alpha=0.7)
    plt.text(z_830 + 0.03, float(phi1[idx]), 'D_C=830 Mpc at z~' + str(round(z_830, 3)), fontsize=9)

    plt.xlabel('Redshift z')
    plt.ylabel('phi(z) (normalized)')
    plt.title('Teardrop Geometry: Profile Comparison')
    plt.grid(True, alpha=0.25)
    plt.legend(frameon=False)
    plt.tight_layout()

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220, bbox_inches='tight')
    fig.savefig(out_pdf, bbox_inches='tight')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rs1', type=float, default=830.0)
    ap.add_argument('--rs2', type=float, default=2000.0)
    ap.add_argument('--outdir', type=str, default='figures')
    args = ap.parse_args()

    outdir = Path(args.outdir)
    out_png = outdir / 'figureA1_teardrop_geometry.png'
    out_pdf = outdir / 'figureA1_teardrop_geometry.pdf'

    make_figure(args.rs1, args.rs2, out_png, out_pdf)


if __name__ == '__main__':
    main()

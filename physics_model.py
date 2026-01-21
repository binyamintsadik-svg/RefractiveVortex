"""physics_model.py

Minimal reproducibility module for the refractive-vortex cosmology calculations used in this notebook.

This is written to be readable and easy to lift into a paper repo. It implements:
- Flat LCDM E(z)
- Comoving distance integrals
- Refractive index modifier eta(z, nu) = 1 + xi(nu) * |phi(x)|, with phi(x) = -ln(1+x)/x and x=r/r_s
- SN distance modulus prediction, plus intrinsic refractive magnitude term gamma*(phi/phi_anchor - 1)
- CMB acoustic-scale proxy l_A = pi * D_M(z*) / r_s(Planck) with a single constraint

Notes:
- r_s in this model is the refractive scale (fixed at 2000 in the analysis), not the baryon sound horizon.
- Frequency chromatic scaling uses xi(nu) = xi0 * nu^alpha, with nu_CMB as a large fiducial.
"""

import numpy as np


def cumtrapz(y_vals, x_vals):
    y_vals = np.asarray(y_vals, dtype=float)
    x_vals = np.asarray(x_vals, dtype=float)
    out = np.zeros_like(y_vals)
    dx = np.diff(x_vals)
    out[1:] = np.cumsum(0.5 * (y_vals[1:] + y_vals[:-1]) * dx)
    return out


def Ez_lcdm(z_vals, Om=0.3):
    z_vals = np.asarray(z_vals, dtype=float)
    Ol = 1.0 - float(Om)
    return np.sqrt(float(Om) * (1.0 + z_vals) ** 3 + Ol)


def phi_from_x(x_vals):
    x_vals = np.asarray(x_vals, dtype=float)
    return -np.log1p(x_vals) / (x_vals + 1e-12)


def build_integrals(z_grid, Om=0.3):
    z_grid = np.asarray(z_grid, dtype=float)
    Ez = Ez_lcdm(z_grid, Om=Om)
    I_Einv = cumtrapz(1.0 / Ez, z_grid)
    return Ez, I_Einv


def phi_arrays_from_integral(I_sn, I_anchor, I_grid, H_base, r_s, c_light):
    coef = (float(c_light) / float(H_base)) / float(r_s)
    phi_sn = phi_from_x(coef * np.asarray(I_sn, dtype=float))
    phi_anchor = float(phi_from_x(coef * float(I_anchor)))
    phi_grid = phi_from_x(coef * np.asarray(I_grid, dtype=float))
    return phi_sn, phi_anchor, phi_grid


def mu_model_sn(z_sn, z_grid, Ez_grid, I_grid, I_sn, H_base, r_s, c_light, xi0, gamma, z_anchor=0.05):
    z_sn = np.asarray(z_sn, dtype=float)

    # Phi arrays
    I_anchor = float(np.interp(float(z_anchor), z_grid, I_grid))
    phi_sn, phi_anchor, phi_grid = phi_arrays_from_integral(I_sn=I_sn, I_anchor=I_anchor, I_grid=I_grid,
                                                            H_base=H_base, r_s=r_s, c_light=c_light)

    # Refractive eta for optical band (nu factor = 1)
    eta_grid = 1.0 + float(xi0) * np.abs(phi_grid)
    J = cumtrapz(1.0 / (Ez_grid * eta_grid), z_grid)
    J_at = np.interp(z_sn, z_grid, J)

    Dc = (float(c_light) / float(H_base)) * J_at
    Dl = (1.0 + z_sn) * Dc

    mu_geom = 5.0 * np.log10(Dl) + 25.0
    mu_intr = float(gamma) * (phi_sn / (phi_anchor + 1e-15) - 1.0)

    return mu_geom + mu_intr


def lA_model(z_grid, Ez_grid, I_grid, H_base, r_s, c_light, xi0, alpha, nu_cmb, r_s_planck):
    # phi grid and eta for CMB frequency
    I_anchor = float(np.interp(0.05, z_grid, I_grid))
    _, _, phi_grid = phi_arrays_from_integral(I_sn=np.array([I_anchor]), I_anchor=I_anchor, I_grid=I_grid,
                                              H_base=H_base, r_s=r_s, c_light=c_light)

    xi_cmb = float(xi0) * float(nu_cmb) ** float(alpha)
    eta_grid = 1.0 + xi_cmb * np.abs(phi_grid)
    J = cumtrapz(1.0 / (Ez_grid * eta_grid), z_grid)
    Dm = (float(c_light) / float(H_base)) * float(J[-1])
    return float(np.pi * Dm / float(r_s_planck))

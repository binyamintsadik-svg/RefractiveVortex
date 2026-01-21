"""main_analysis.py

Reproducibility entry-point:
- Loads Pantheon+SH0ES-like SN data (expects a CSV with columns: z, mu, mu_err)
- Computes SN+CMB chi2
- Runs a grid scan over (H_base, alpha) by profiling over (xi0, gamma)
- Runs constrained Monte Carlo with alpha fixed = 1/3 to estimate H_base uncertainty

This is intended as a clean, minimal re-run script. For the exact dataset used in the notebook,
export your SN table into the expected CSV format.
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from tqdm import tqdm

import physics_model as pm

# ---------- User-set inputs ----------
SN_CSV = 'sn_data.csv'   # must contain: z, mu, mu_err

# Cosmology settings
OM = 0.3
C_LIGHT = 299792.458

# Refractive settings
R_S_REF = 2000.0
NU_CMB = 1.0e11

# CMB constraint
R_S_PLANCK = 147.05
L_A_PLANCK = 301.63
SIGMA_LA = 0.18

# Grid scan settings
H_VALS = np.linspace(66.6, 68.4, 41)
ALPHA_VALS = np.linspace(0.22, 0.50, 57)

# Monte Carlo settings
MC_N = 300
ALPHA_FIXED = 1.0 / 3.0


def profile_xi0_gamma(z_sn, mu_obs, mu_err, z_grid, Ez_grid, I_grid, I_sn, H_base, alpha_val, xi0_bounds=(0.0, 10.0), gamma_bounds=(-5.0, 5.0)):
    def chi2_inner(x_pair):
        xi0 = float(x_pair[0])
        gamma = float(x_pair[1])

        mu_pred = pm.mu_model_sn(z_sn=z_sn, z_grid=z_grid, Ez_grid=Ez_grid, I_grid=I_grid, I_sn=I_sn,
                                H_base=H_base, r_s=R_S_REF, c_light=C_LIGHT, xi0=xi0, gamma=gamma)
        chi2_sn = float(np.sum(((mu_obs - mu_pred) / mu_err) ** 2))

        lA = pm.lA_model(z_grid=z_grid, Ez_grid=Ez_grid, I_grid=I_grid, H_base=H_base, r_s=R_S_REF, c_light=C_LIGHT,
                         xi0=xi0, alpha=alpha_val, nu_cmb=NU_CMB, r_s_planck=R_S_PLANCK)
        chi2_cmb = float(((lA - L_A_PLANCK) / SIGMA_LA) ** 2)

        return chi2_sn + chi2_cmb

    res = minimize(chi2_inner, x0=np.array([0.10, -0.2]), method='L-BFGS-B',
                   bounds=[xi0_bounds, gamma_bounds], options={'maxiter': 250})
    return float(res.fun), float(res.x[0]), float(res.x[1])


def main():
    sn_df = pd.read_csv(SN_CSV)
    z_sn = sn_df['z'].to_numpy(dtype=float)
    mu_obs = sn_df['mu'].to_numpy(dtype=float)
    mu_err = sn_df['mu_err'].to_numpy(dtype=float)

    # grids
    z_grid = np.unique(np.concatenate([np.linspace(0.0, 2.5, 2000), np.logspace(-4, np.log10(1100.0), 8000)]))
    z_grid.sort()
    Ez_grid, I_grid = pm.build_integrals(z_grid, Om=OM)
    I_sn = np.interp(z_sn, z_grid, I_grid)

    # grid scan
    chi2_grid = np.zeros((len(ALPHA_VALS), len(H_VALS)))
    xi0_grid = np.zeros_like(chi2_grid)
    gamma_grid = np.zeros_like(chi2_grid)

    for ia, a in enumerate(tqdm(ALPHA_VALS)):
        for ih, H in enumerate(H_VALS):
            c2, xi0_b, g_b = profile_xi0_gamma(z_sn, mu_obs, mu_err, z_grid, Ez_grid, I_grid, I_sn, H, a)
            chi2_grid[ia, ih] = c2
            xi0_grid[ia, ih] = xi0_b
            gamma_grid[ia, ih] = g_b

    # Constrained MC
    rng = np.random.default_rng(24680)
    H_samps = np.zeros(MC_N)

    for i in tqdm(range(MC_N)):
        mu_pert = mu_obs + rng.normal(0.0, mu_err)

        def obj_H(H_arr):
            H = float(H_arr[0])
            c2, _, _ = profile_xi0_gamma(z_sn, mu_pert, mu_err, z_grid, Ez_grid, I_grid, I_sn, H, ALPHA_FIXED)
            return c2

        resH = minimize(obj_H, x0=np.array([67.5]), method='L-BFGS-B', bounds=[(63.0, 72.0)], options={'maxiter': 150})
        H_samps[i] = float(resH.x[0])

    print('H_base mean')
    print(float(np.mean(H_samps)))
    print('H_base std')
    print(float(np.std(H_samps, ddof=1)))

    np.savez('analysis_outputs.npz', H_vals=H_VALS, alpha_vals=ALPHA_VALS, chi2_grid=chi2_grid,
             xi0_grid=xi0_grid, gamma_grid=gamma_grid, H_mc=H_samps)


if __name__ == '__main__':
    main()

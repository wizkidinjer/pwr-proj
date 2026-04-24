# setup_checker.py
"""
Sweep design variables, find min V_assy config that delivers P = 500 MW.

Tuning knobs at top: edit the linspaces, run the file.
"""
import numpy as np
import matplotlib.pyplot as plt

import config as c
from core import Core
from plotter import plot_core_xy


# =========================================================================
# SWEEP RANGES -- edit these
# =========================================================================
D_ROD_RANGE      = np.linspace(8.0e-3, 11.0e-3, 5)     # m, rod outer dia
ASSY_PITCH_RANGE = np.linspace(0.15,    0.25,    5)    # m, assembly edge
N_SIDE_RANGE     = np.array([15, 17, 19, 21])          # pins per assy edge
N_ASSY_RANGE     = [5, 9, 21, 25]                      # layouts available

D_ROD_RANGE = [10.25e-3]
ASSY_PITCH_RANGE=[.2]
N_SIDE_RANGE = [19]
N_ASSY_RANGE=[25]
# =========================================================================
# OPTIONAL FILTERS (set None to disable)
# =========================================================================
P_OVER_D_MIN = 1.15     # skip configs tighter than this (coolant gap)
H_MAX        = None     # m, cap on active length (e.g. 4.5 for realism)


def run_sweep():
    best = None   # (V_assy, D_rod, p_assy, N_side, N_assy, H, P_check)
    n_total = 0
    n_valid = 0

    for N_assy in N_ASSY_RANGE:
        for N_side in N_SIDE_RANGE:
            for p_assy in ASSY_PITCH_RANGE:
                rod_pitch = p_assy / N_side
                for D_rod in D_ROD_RANGE:
                    n_total += 1

                    # --- validity filters ---
                    if D_rod >= rod_pitch:
                        continue
                    if P_OVER_D_MIN is not None:
                        if rod_pitch / D_rod < P_OVER_D_MIN:
                            continue

                    # --- build core, solve H ---
                    try:
                        core = Core(
                            D_rod=D_rod,
                            assembly_pitch=p_assy,
                            N_assemblies=N_assy,
                            N_side=int(N_side),
                            clad_thick=c.CLAD_THICK,
                            H_active=1.0,   # placeholder
                        )
                    except ValueError:
                        continue

                    H = core.solve_H_for_power(
                        P_target=c.P_THERMAL,
                        phi_max=c.PHI_MAX,
                        N_ff=c.N_FF,
                        sigma_f=c.SIGMA_F,
                        G_f=c.G_F_J,
                    )

                    if H_MAX is not None and H > H_MAX:
                        continue

                    # --- solve mass flow per assembly ---
                    assy = core.assemblies[0]
                    assy.MAX_P = c.DP_CORE_MAX
                    assy.P_DIFF = c.DP_CORE_MAX
                    mdot_assy, *_ = assy.iterate_m_dot(1000,verbose=True)  # initial guess = 1 kg/s
                    mdot_core = mdot_assy * N_assy

                    V_assy = core.V_assemblies
                    n_valid += 1

                    if best is None or V_assy < best[0]:
                        P_check = core.total_power(
                            c.PHI_MAX, c.N_FF, c.SIGMA_F, c.G_F_J
                        )
                        best = (V_assy, D_rod, p_assy, int(N_side),
                                N_assy, H, P_check, core,
                                mdot_assy, mdot_core)

    print(f"\nSwept {n_total} configs, {n_valid} valid.\n")

    if best is None:
        print("No valid config. Loosen filters.")
        return None

    V, D_rod, p_assy, N_side, N_assy, H, P_check, core, mdot_assy, mdot_core = best
    dP_MW = (P_check - c.P_THERMAL) / 1e6

    print("=" * 60)
    print("MIN VOLUME CONFIG")
    print("=" * 60)
    print(f"  V_assemblies    = {V:.4f} m^3   (OPTIMIZATION TARGET)")
    print(f"  P_total (check) = {P_check/1e6:.4f} MW   (target 500.0000)")
    print(f"  |P - 500 MW|    = {abs(dP_MW):.2e} MW   (tol 0.1)")
    print("-" * 60)
    print(f"  D_rod           = {D_rod*1e3:.3f} mm")
    print(f"  assembly_pitch  = {p_assy*1e3:.2f} mm")
    print(f"  N_side          = {N_side}  ({N_side**2} rods/assy)")
    print(f"  N_assemblies    = {N_assy}  ({N_assy * N_side**2} rods total)")
    print(f"  H_active        = {H:.4f} m")
    print("-" * 60)
    print(f"  rod_pitch       = {p_assy/N_side*1e3:.3f} mm")
    print(f"  P/D             = {(p_assy/N_side)/D_rod:.3f}")
    print(f"  rod_spacing     = {(p_assy/N_side - D_rod)*1e3:.3f} mm")
    print(f"  R_core          = {core.R_core*1e3:.2f} mm")
    print(f"  D_core          = {2*core.R_core:.4f} m")
    print(f"  V_core_cyl      = {core.V_core_cyl:.4f} m^3  (reference)")
    print("-" * 60)
    print(f"  mdot per assy   = {mdot_assy:.4f} kg/s")
    print(f"  mdot total core = {mdot_core:.4f} kg/s")
    print("=" * 60)

    return core


if __name__ == "__main__":
    core = run_sweep()
    # if core is not None:
    #     plot_core_xy(core, savepath="min_vol_config.png")
    #     plt.show()
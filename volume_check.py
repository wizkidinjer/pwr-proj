# volume_check.py
"""
Sweep (D_rod, P/D, N_side, N_assy) -> find min V_total config.
Pure geometry. H is forced by physics (500 MW). No fluid math here.

Output: prints min config, saves best_config.json for main.py.
"""
import json
import numpy as np
from tqdm import tqdm
import config as c
from geometry import Core, LAYOUTS

# ============ SWEEP RANGES ============
D_ROD_RANGE  = np.linspace(1.0e-3, 10.0e-3, 3)
P_OVER_D_RNG = np.linspace(1.0,   1.4,  30)
N_SIDE_RANGE = [17]
N_ASSY_RANGE = sorted(LAYOUTS)            # all 9 layouts

# Optional: skip configs above this rod length for realism (set None to disable)
H_MAX = None


def run_sweep():
    best = None
    n_total, n_valid = 0, 0
    for N_assy in tqdm(N_ASSY_RANGE):
        for N_side in N_SIDE_RANGE:
            for P_over_D in P_OVER_D_RNG:
                for D_rod in D_ROD_RANGE:
                    n_total += 1
                    if P_over_D <= 1.0:
                        continue
                    try:
                        core = Core(D_rod, P_over_D, N_side, N_assy, H=1.0)
                    except ValueError:
                        continue
                    H = core.solve_H()
                    if H_MAX is not None and H > H_MAX:
                        continue
                    n_valid += 1
                    V = core.V_total
                    if best is None or V < best["V_total"]:
                        best = {
                            "D_rod":     float(D_rod),
                            "P_over_D":  float(P_over_D),
                            "N_side":    int(N_side),
                            "N_assy":    int(N_assy),
                            "p_assy":    float(core.p_assy),
                            "rod_pitch": float(core.rod_pitch),
                            "H":         float(H),
                            "R_core":    float(core.R_core),
                            "V_total":   float(V),
                            "rods_per_assy": int(N_side ** 2),
                            "rods_total":    int(core.n_pins_total),
                        }
    return best, n_total, n_valid


def main():
    best, n_total, n_valid = run_sweep()
    print(f"\nSwept {n_total} configs, {n_valid} valid.\n")
    if best is None:
        print("No valid config. Loosen filters."); return

    print("=" * 56)
    print("MIN-VOLUME CONFIG")
    print("=" * 56)
    print(f"  V_total        = {best['V_total']:.4f} m^3   (target)")
    print(f"  H_active       = {best['H']:.4f} m")
    print("-" * 56)
    print(f"  D_rod          = {best['D_rod']*1e3:.3f} mm")
    print(f"  P/D            = {best['P_over_D']:.3f}")
    print(f"  rod_pitch      = {best['rod_pitch']*1e3:.3f} mm")
    print(f"  N_side         = {best['N_side']}  ({best['rods_per_assy']} rods/assy)")
    print(f"  p_assy         = {best['p_assy']*1e3:.2f} mm")
    print(f"  N_assy         = {best['N_assy']}  ({best['rods_total']} rods total)")
    print(f"  R_core         = {best['R_core']*1e3:.2f} mm")
    print(f"  D_core         = {2*best['R_core']:.4f} m")
    print(f"  rod_spacing    = {(best['rod_pitch']-best['D_rod'])*1e3:.3f} mm")
    print(f"  filled_frac    = {( (best['p_assy']**2*9) / (best['R_core']**2*np.pi) ) :.3f}")
    print("=" * 56)

    with open("best_config.json", "w") as f:
        json.dump(best, f, indent=2)
    print(f"\nSaved to best_config.json")

    # rebuild core for plotting
    from geometry import Core
    from plotter import plot_core_xy, plot_core_3d
    import matplotlib.pyplot as plt

    core = Core(best["D_rod"], best["P_over_D"],
                best["N_side"], best["N_assy"])
    core.solve_H()

    plot_core_xy(core, savepath="best_xy.png")
    plot_core_3d(core, show_rods=False, savepath="best_3d.png")
    plt.show()


if __name__ == "__main__":
    main()

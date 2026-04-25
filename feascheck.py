# feasible_check.py
"""
Sweep (D_rod, P/D, N_side, N_assy). Run thermal solve on each.
Track min V_total among configs that PASS both:
  (1) T_out + T_SAT_MARGIN <= T_sat(P_sys)
  (2) T_0_max <= T_FUEL_MAX

Saves winner to best_config.json (same format as volume_check.py).
"""
import json
import time
import numpy as np
from tqdm import tqdm
import config as c
import water_props as wp
from geometry import Core, LAYOUTS
from main import solve_flow_and_temps, axial_temps

# ============ SWEEP RANGES ============
D_ROD_RANGE  = np.linspace(1.0e-3, 11.0e-3, 20)
P_OVER_D_RNG = np.linspace(1.0,   1.50,    20)
N_SIDE_RANGE = [15, 17, 19, 21]
N_ASSY_RANGE = sorted(LAYOUTS)


def evaluate(core):
    """Return (passes, V_total, T_out, T_0_max, margin) or (False, ...) if solver fail."""
    try:
        P_rod_c = core.center_pin().power(core.R_core)
        flow  = solve_flow_and_temps(core, P_rod_c)
        therm = axial_temps(core, flow, P_rod_c)
        T_sat = wp.T_sat(c.P_SYSTEM)
        margin = T_sat - flow["T_out"]
        ok = (margin >= c.T_SAT_MARGIN) and (therm["T_0_max"] <= c.T_FUEL_MAX)
        return ok, core.V_total, flow["T_out"], therm["T_0_max"], margin
    except Exception:
        return False, None, None, None, None


def run_sweep():
    best, best_record = None, None
    n_total = n_built = n_pass = 0
    t0 = time.time()
    for N_assy in tqdm(N_ASSY_RANGE):
        for N_side in N_SIDE_RANGE:
            for P_over_D in P_OVER_D_RNG:
                for D_rod in D_ROD_RANGE:
                    n_total += 1
                    if P_over_D <= 1.0:
                        continue
                    try:
                        core = Core(D_rod, P_over_D, N_side, N_assy)
                    except ValueError:
                        continue
                    core.solve_H()
                    n_built += 1
                    ok, V, T_out, T0_max, margin = evaluate(core)
                    if not ok:
                        continue
                    n_pass += 1
                    if best is None or V < best:
                        best = V
                        best_record = {
                            "D_rod":     float(D_rod),
                            "P_over_D":  float(P_over_D),
                            "N_side":    int(N_side),
                            "N_assy":    int(N_assy),
                            "p_assy":    float(core.p_assy),
                            "rod_pitch": float(core.rod_pitch),
                            "H":         float(core.H),
                            "R_core":    float(core.R_core),
                            "V_total":   float(V),
                            "T_out":     float(T_out),
                            "T_0_max":   float(T0_max),
                            "margin":    float(margin),
                            "rods_per_assy": int(N_side**2),
                            "rods_total":    int(core.n_pins_total),
                        }
    return best_record, n_total, n_built, n_pass, time.time()-t0


def main():
    rec, n_total, n_built, n_pass, dt = run_sweep()
    print(f"\nSwept {n_total}, built {n_built}, passing {n_pass}.  ({dt:.1f}s)\n")
    if rec is None:
        print("No config passes both constraints. Loosen sweep or accept failure.")
        return

    print("=" * 60)
    print("MIN-VOLUME FEASIBLE CONFIG")
    print("=" * 60)
    print(f"  V_total       = {rec['V_total']:.4f} m^3")
    print(f"  H_active      = {rec['H']:.4f} m")
    print("-" * 60)
    print(f"  D_rod         = {rec['D_rod']*1e3:.3f} mm")
    print(f"  P/D           = {rec['P_over_D']:.3f}")
    print(f"  N_side        = {rec['N_side']}  ({rec['rods_per_assy']} rods/assy)")
    print(f"  N_assy        = {rec['N_assy']}  ({rec['rods_total']} rods total)")
    print(f"  p_assy        = {rec['p_assy']*1e3:.2f} mm")
    print(f"  R_core        = {rec['R_core']*1e3:.2f} mm")
    print("-" * 60)
    print(f"  T_out         = {rec['T_out']:.2f} C")
    print(f"  margin        = {rec['margin']:.2f} K   (req >= {c.T_SAT_MARGIN})")
    print(f"  T_0_max       = {rec['T_0_max']:.2f} C  (limit {c.T_FUEL_MAX})")
    print("=" * 60)

    with open("best_config.json", "w") as f:
        json.dump(rec, f, indent=2)
    print("\nSaved to best_config.json (run main.py for full report + plot)")


if __name__ == "__main__":
    main()
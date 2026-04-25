# main.py
"""
Read best_config.json, run full thermal-hydraulic solve:
  - iterate (mdot, T_mean) jointly using IAPWS water props
  - axial T_c(z), T_0(z) profiles, find hot spot
  - check T_sat margin and T_fuel_max constraints
"""
import json
import numpy as np
import config as c
import water_props as wp
from geometry import Core


# ----------- friction (Cheng-Todreas, square lattice turbulent) -----------
def friction_factor(Re, P_over_D):
    C = c.CT_A + c.CT_B1*(P_over_D-1) + c.CT_B2*(P_over_D-1)**2
    return C / Re**0.18


# ----------- coupled (mdot, T_mean) fixed-point -----------
def solve_flow_and_temps(core, P_rod_center,
                         dp=c.DP_CORE_MAX, T_in=c.T_IN, P_sys=c.P_SYSTEM,
                         tol=1e-4, max_iter=50):
    """
    Returns dict with mdot_assy, mdot_rod, T_out, T_mean, props (at T_mean).
    Outer loop: T_mean -> props -> mdot, T_out -> new T_mean.
    """
    p_assy, P_over_D = core.p_assy, core.P_over_D
    D, H, N_pins = core.D_rod, core.H, core.N_side**2
    pin_pitch = core.rod_pitch
    A_c = pin_pitch**2 - np.pi*D**2/4               # subchannel area per rod
    A_assy = A_c * N_pins                            # total flow area per assembly
    D_h = 4 * A_c / (np.pi * D)                      # rod-channel hydraulic dia

    T_mean = T_in
    for it in range(max_iter):
        pr = wp.props(T_mean, P_sys)
        # mdot per assy: solve dp = f * H/D_h * rho u^2 / 2,  u = mdot_assy/(rho A_assy)
        # closed form vs Re given Cheng-Todreas: iterate Re because f depends on Re
        mdot = 1.0
        for _ in range(50):
            u  = mdot / (pr["rho"] * A_assy)
            Re = pr["rho"] * u * D_h / pr["mu"]
            f  = friction_factor(Re, P_over_D)
            mdot_new = np.sqrt(dp * 2 * pr["rho"] * A_assy**2 / (f * H / D_h))
            if abs(mdot_new - mdot) < 1e-6 * mdot:
                mdot = mdot_new; break
            mdot = mdot_new
        mdot_rod = mdot / N_pins
        T_out = T_in + P_rod_center / (mdot_rod * pr["cp"])
        T_mean_new = 0.5 * (T_in + T_out)
        if abs(T_mean_new - T_mean) < tol:
            T_mean = T_mean_new; break
        T_mean = T_mean_new

    return {
        "props": pr, "mdot_assy": mdot, "mdot_rod": mdot_rod,
        "T_in": T_in, "T_out": T_out, "T_mean": T_mean,
        "Re": Re, "u": u, "D_h": D_h, "A_c": A_c, "iters": it+1,
    }


# ----------- axial temperature profile, hot-spot search -----------
def axial_temps(core, flow, P_rod_center, n_z=2000):
    """Vectorized axial profiles for the center rod. Returns z, T_c, T_co, T_ci, T_fo, T_0."""
    H, D = core.H, core.D_rod
    D_fuel = D - 2*c.CLAD_THICK
    r_co, r_ci = D/2, D/2 - c.CLAD_THICK
    pr   = flow["props"]
    mdot = flow["mdot_rod"]

    # Dittus-Boelter h (heating: Pr^0.4)
    Nu = 0.023 * flow["Re"]**0.8 * pr["pr"]**0.4
    h  = Nu * pr["k"] / flow["D_h"]

    # resistances per unit length [m K / W]
    R_conv = 1 / (np.pi * D     * h)
    R_clad = np.log(r_co/r_ci) / (2 * np.pi * c.K_CLAD)
    R_gap  = 1 / (np.pi * D_fuel * c.H_GAP)
    R_fuel = 1 / (4 * np.pi * c.K_FUEL)
    R_tot  = R_conv + R_clad + R_gap + R_fuel

    # peak linear power for center pin (J0 = 1)
    A_fuel = np.pi * (D_fuel/2)**2
    qp0    = c.N_FF * c.SIGMA_F * c.G_F_J * c.PHI_MAX * A_fuel   # W/m

    z   = np.linspace(-H/2, H/2, n_z)
    qp  = qp0 * np.cos(np.pi * z / H)
    # T_c from energy balance (integral of cos = (H/pi)(sin+1))
    T_c  = c.T_IN + qp0 * H / (mdot * pr["cp"] * np.pi) * (np.sin(np.pi*z/H) + 1)
    T_co = T_c  + qp * R_conv
    T_ci = T_co + qp * R_clad
    T_fo = T_ci + qp * R_gap
    T_0  = T_fo + qp * R_fuel

    idx_hot = int(np.argmax(T_0))
    return {
        "z": z, "T_c": T_c, "T_co": T_co, "T_ci": T_ci, "T_fo": T_fo, "T_0": T_0,
        "idx_hot": idx_hot, "z_hot": float(z[idx_hot]),
        "T_0_max": float(T_0[idx_hot]),
        "T_c_at_hot": float(T_c[idx_hot]),
        "h_conv": float(h), "qp0": float(qp0),
        "R": {"conv": R_conv, "clad": R_clad, "gap": R_gap, "fuel": R_fuel, "tot": R_tot},
    }


# ----------- main -----------
def main(plot=False):
    with open("best_config.json") as f:
        cfg = json.load(f)
    core = Core(cfg["D_rod"], cfg["P_over_D"], cfg["N_side"], cfg["N_assy"])
    core.solve_H()
    P_check = core.total_power()
    P_rod_c = core.center_pin().power(core.R_core)

    flow = solve_flow_and_temps(core, P_rod_c)
    therm = axial_temps(core, flow, P_rod_c)
    T_sat = wp.T_sat(c.P_SYSTEM)

    margin   = T_sat - flow["T_out"]
    pass_sat = margin >= c.T_SAT_MARGIN
    pass_fuel = therm["T_0_max"] <= c.T_FUEL_MAX

    print("=" * 60)
    print("CORE (from best_config.json)")
    print("=" * 60)
    print(f"  D_rod         {core.D_rod*1e3:.3f} mm     "
          f"P/D = {core.P_over_D:.3f}")
    print(f"  N_side        {core.N_side}        "
          f"N_assy = {core.N_assy}")
    print(f"  p_assy        {core.p_assy*1e3:.2f} mm")
    print(f"  H_active      {core.H:.4f} m   "
          f"V_total = {core.V_total:.4f} m^3")
    print(f"  R_core        {core.R_core*1e3:.2f} mm   "
          f"P_total(check) = {P_check/1e6:.4f} MW")

    print("\n" + "=" * 60)
    print("FLOW (coupled mdot + T_mean iteration)")
    print("=" * 60)
    pr = flow["props"]
    print(f"  iters         {flow['iters']}")
    print(f"  T_mean        {flow['T_mean']:.2f} C   props at T_mean:")
    print(f"     rho={pr['rho']:.1f}  mu={pr['mu']:.3e}  k={pr['k']:.4f}  "
          f"cp={pr['cp']:.0f}  Pr={pr['pr']:.3f}")
    print(f"  mdot per assy {flow['mdot_assy']:.3f} kg/s")
    print(f"  mdot per rod  {flow['mdot_rod']:.4f} kg/s")
    print(f"  Re            {flow['Re']:.3e}")
    print(f"  T_in / T_out  {flow['T_in']:.2f}  /  {flow['T_out']:.2f}  C")

    print("\n" + "=" * 60)
    print("HOT SPOT (center rod)")
    print("=" * 60)
    print(f"  qp0           {therm['qp0']:.1f} W/m  (peak linear power)")
    print(f"  h_conv        {therm['h_conv']:.0f} W/m^2-K")
    print(f"  z_hot         {therm['z_hot']:.4f} m   (z=0 mid, +H/2 outlet)")
    print(f"  T_c at hot    {therm['T_c_at_hot']:.2f} C")
    print(f"  T_0 max       {therm['T_0_max']:.2f} C")

    print("\n" + "=" * 60)
    print("CONSTRAINTS")
    print("=" * 60)
    print(f"  T_sat at {c.P_SYSTEM/1e6:.1f} MPa = {T_sat:.2f} C")
    print(f"  T_out                  = {flow['T_out']:.2f} C")
    print(f"  margin = T_sat - T_out = {margin:.2f} K   (req >= {c.T_SAT_MARGIN}) "
          f"  {'PASS' if pass_sat else 'FAIL'}")
    print(f"  T_0_max                = {therm['T_0_max']:.2f} C  "
          f"(limit {c.T_FUEL_MAX})        {'PASS' if pass_fuel else 'FAIL'}")
    print("=" * 60)

    if plot:
        import matplotlib.pyplot as plt
        z = therm["z"]
        plt.figure(figsize=(8, 5))
        for label in ("T_c", "T_co", "T_ci", "T_fo", "T_0"):
            plt.plot(z, therm[label], label=label)
        plt.axvline(therm["z_hot"], ls="--", c="k", alpha=0.4, label=f"z_hot={therm['z_hot']:.3f}")
        plt.axhline(c.T_FUEL_MAX, ls=":", c="r", label=f"T_fuel_max={c.T_FUEL_MAX}")
        plt.xlabel("z (m)"); plt.ylabel("T (C)"); plt.legend(); plt.grid(alpha=0.3)
        plt.title("Center rod axial temperatures")
        plt.savefig("axial_temps.png", dpi=140, bbox_inches="tight")
        plt.show()


if __name__ == "__main__":
    main(plot=True)

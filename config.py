# config.py
"""Constants only. Tuneable design vars live in volume_check.py sweep ranges."""

# ------ rubric physics ------
P_THERMAL    = 500e6
PHI_MAX      = 5e13 * 1e4
N_FF         = 8e20 * 1e6
SIGMA_F      = 585e-28
G_F_J        = 200.0 * 1.602176634e-13

# ------ PWR operating ------
P_SYSTEM     = 15.5e6
DP_CORE_MAX  = 175e3
T_FUEL_MAX   = 1000.0          # deg C, max centerline
T_IN         = 290.0           # deg C, inlet coolant
T_SAT_MARGIN = 20.0            # K, required (T_sat - T_out)

# ------ materials ------
H_GAP        = 7000.0
K_CLAD       = 20.0
K_FUEL       = 2.0
CLAD_THICK   = 0.65e-3

# ------ Cheng-Todreas friction (square lattice turbulent) ------
CT_A, CT_B1, CT_B2 = 0.1339, 0.09059, -0.09926

J0_FIRST_ZERO = 2.405

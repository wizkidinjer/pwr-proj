# config.py
"""
Central config for PWR core thermal-hydraulic design.
Two blocks:
  1. PHYSICS / RUBRIC CONSTANTS  -- fixed by problem statement, do NOT tune.
  2. DESIGN VARIABLES            -- tune these to hit specs + min core volume.

All lengths in SI (m) unless noted. Rubric values converted where needed.
"""

# =========================================================================
# 1. PHYSICS / RUBRIC CONSTANTS  (locked by problem statement)
# =========================================================================

# --- Core power + neutronics ---
P_THERMAL       = 500e6          # W, core thermal power
PHI_MAX         = 5e13 * 1e4     # neutrons / m^2 / s  (5e13 n/cm^2-s -> m^2)
N_FF            = 8e20 * 1e6     # fissionable nuclei / m^3  (cm^-3 -> m^-3)
SIGMA_F         = 585e-28        # m^2, effective fission xs (585 barn)
G_F_MEV         = 200.0          # MeV per fission
G_F_J           = G_F_MEV * 1.602176634e-13   # J per fission

# Bessel zero used in radial flux shape J0(2.405 r / R)
J0_FIRST_ZERO   = 2.405

# --- PWR operating conditions (rubric table) ---
P_SYSTEM        = 15.5e6         # Pa, nominal core pressure
DP_CORE_MAX     = 175e3          # Pa, max allowable core pressure drop
T_FUEL_MAX      = 1000.0         # deg C, max centerline fuel temp

# --- PWR material / interface properties (rubric table) ---
H_GAP           = 7000.0         # W/m^2-K, fuel-clad gas gap conductance
K_CLAD          = 20.0           # W/m-K, clad thermal conductivity
CLAD_THICK      = 0.65e-3        # m, cladding thickness

# --- Coolant inlet conditions (typical PWR, adjust if spec changes) ---
T_IN            = 290.0          # deg C, core inlet coolant temp
# outlet / mass flow are DERIVED from power balance, not set here.

#water properties
MOD_RHO = 1000
MOD_MU = 0.0010518
MOD_PR = 0.7
MOD_K = 0.564



# =========================================================================
# 2. DESIGN VARIABLES  (tune these)
# =========================================================================

# --- Rod geometry ---
# Pellet OD ~ typical PWR 8.2 mm -> rod OD ~ pellet + 2*clad ~ 9.5 mm.
# Here D_OUTER = fuel-pellet + cladding outer diameter (matches Pin class).
D_ROD_OUTER     = 9.5e-3         # m, rod outer diameter (fuel + clad)

# --- Assembly geometry ---
# Square lattice. pin_pitch = center-to-center spacing.
# Typical PWR P/D ~ 1.30-1.35. Spacing = pitch - D_outer (derived, not stored).
PIN_PITCH       = 12.6e-3        # m, pin-to-pin pitch in assembly
N_SIDE          = 17             # pins per edge (17x17 = 289 rods/assembly)

# --- Core layout ---
# Number of assemblies to place in core. Rubric-style options:
#   5  -> plus (+) pattern
#   9  -> 3x3 square
#   21 -> 3x3 + edge extensions (diamond)
#   25 -> 5x5 square
#   etc.
N_ASSEMBLIES    = 25

# --- Active fuel column height (initial guess, iterated vs volume target) ---
H_ACTIVE        = 3.66           # m, active fuel length (~ref PWR 12 ft)

# --- Extrapolated dimensions for flux shape ---
# Rubric says ignore extrapolation lengths, so R_ex = R_core, H_ex = H_active.
# Kept as separate knobs in case you want to add extrapolation later.
def R_EXTRAPOLATED(R_core):
    return R_core     # no extrapolation per rubric
def H_EXTRAPOLATED(H_active):
    return H_active   # no extrapolation per rubric


# =========================================================================
# 3. DERIVED HELPERS (read-only, computed from above)
# =========================================================================

def rod_spacing():
    """Gap between adjacent rod surfaces (m). pitch - D_outer."""
    return PIN_PITCH - D_ROD_OUTER

def p_over_d():
    """Pitch-to-diameter ratio. Typical PWR 1.30-1.35."""
    return PIN_PITCH / D_ROD_OUTER

def assembly_outer():
    """Assembly outer edge length (m). Square, = N_side * pitch."""
    return N_SIDE * PIN_PITCH

def n_rods_per_assembly():
    return N_SIDE ** 2

def n_rods_total():
    return N_ASSEMBLIES * n_rods_per_assembly()
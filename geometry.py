# geometry.py
"""Pin, Assembly, Core in one file. Vectorized hot paths."""
import numpy as np
from scipy.special import j0
import config as c


# -------------------- layouts --------------------
def _gauss_circle(N_target):
    for R in np.linspace(0, 10, 8000):
        rng = int(np.ceil(R)) + 1
        cells = [(i, j) for i in range(-rng, rng + 1)
                 for j in range(-rng, rng + 1) if np.hypot(i, j) <= R]
        if len(cells) == N_target:
            return sorted(cells)
    raise ValueError(f"No Gauss-circle layout of size {N_target}.")


LAYOUTS = {
    5:  [(0,0),(1,0),(-1,0),(0,1),(0,-1)],
    9:  [(i,j) for i in (-1,0,1) for j in (-1,0,1)],
    21: [(i,j) for i in (-1,0,1) for j in (-1,0,1)]
        + [(2,0),(-2,0),(0,2),(0,-2)]
        + [(2,1),(2,-1),(-2,1),(-2,-1),(1,2),(-1,2),(1,-2),(-1,-2)],
    25: [(i,j) for i in range(-2,3) for j in range(-2,3)],
    37: _gauss_circle(37), 57: _gauss_circle(57), 69: _gauss_circle(69),
    89: _gauss_circle(89), 97: _gauss_circle(97),
}


# -------------------- Pin --------------------
class Pin:
    __slots__ = ("D_outer", "clad_thick", "D_fuel", "H", "x", "y")

    def __init__(self, D_outer, clad_thick, H, xy):
        self.D_outer, self.clad_thick = D_outer, clad_thick
        self.D_fuel = D_outer - 2 * clad_thick
        if self.D_fuel <= 0:
            raise ValueError(f"Clad too thick for D_outer={D_outer}")
        self.H, (self.x, self.y) = H, xy

    @property
    def r(self):       return np.hypot(self.x, self.y)
    @property
    def A_fuel(self):  return np.pi * (self.D_fuel / 2) ** 2
    @property
    def A_clad(self):  return np.pi * ((self.D_outer/2)**2 - (self.D_fuel/2)**2)
    @property
    def V_fuel(self):  return self.A_fuel * self.H
    @property
    def V_clad(self):  return self.A_clad * self.H

    def power(self, R_core, phi_max=c.PHI_MAX, N_ff=c.N_FF,
              sigma_f=c.SIGMA_F, G_f=c.G_F_J):
        """Rod power [W] under thin-rod assumption."""
        return (N_ff * sigma_f * G_f * phi_max
                * j0(c.J0_FIRST_ZERO * self.r / R_core)
                * self.A_fuel * 2 * self.H / np.pi)


# -------------------- Assembly --------------------
class Assembly:
    def __init__(self, N_side, p_assy, center_xy, D_rod, clad_thick, H):
        self.N_side, self.p_assy = N_side, p_assy
        self.pin_pitch = p_assy / N_side
        if D_rod >= self.pin_pitch:
            raise ValueError(f"D={D_rod*1e3:.2f}mm >= pitch={self.pin_pitch*1e3:.2f}mm")
        self.cx, self.cy = center_xy
        self.D_rod, self.H = D_rod, H
        # build pins
        x0 = self.cx - p_assy/2 + self.pin_pitch/2
        y0 = self.cy - p_assy/2 + self.pin_pitch/2
        self.pins = [Pin(D_rod, clad_thick, H,
                         (x0 + i*self.pin_pitch, y0 + j*self.pin_pitch))
                     for i in range(N_side) for j in range(N_side)]

    @property
    def n_pins(self):       return self.N_side ** 2
    @property
    def rod_spacing(self):  return self.pin_pitch - self.D_rod
    @property
    def bbox(self):
        h = self.p_assy / 2
        return (self.cx-h, self.cy-h, self.cx+h, self.cy+h)
    @property
    def V_coolant(self):
        return self.p_assy**2 * self.H - sum(p.V_fuel + p.V_clad for p in self.pins)


# -------------------- Core --------------------
class Core:
    def __init__(self, D_rod, P_over_D, N_side, N_assy,
                 clad_thick=c.CLAD_THICK, H=1.0):
        if N_assy not in LAYOUTS:
            raise ValueError(f"N_assy={N_assy} not in {sorted(LAYOUTS)}")
        if P_over_D <= 1.0:
            raise ValueError(f"P/D={P_over_D} must be > 1 (rods overlap)")
        self.D_rod, self.P_over_D, self.N_side, self.N_assy = D_rod, P_over_D, N_side, N_assy
        self.rod_pitch = D_rod * P_over_D
        self.p_assy    = N_side * self.rod_pitch
        self.H         = H

        # build assemblies
        offsets = LAYOUTS[N_assy]
        self.assemblies = [
            Assembly(N_side, self.p_assy, (i*self.p_assy, j*self.p_assy),
                     D_rod, clad_thick, H)
            for (i, j) in offsets
        ]

        # vectorized pin positions (for fast J0 sums)
        self._x = np.array([p.x for a in self.assemblies for p in a.pins])
        self._y = np.array([p.y for a in self.assemblies for p in a.pins])
        self._r = np.hypot(self._x, self._y)

        # R_core = max corner distance over all assemblies
        corners = np.array([
            (cx, cy)
            for a in self.assemblies
            for cx in (a.bbox[0], a.bbox[2])
            for cy in (a.bbox[1], a.bbox[3])
        ])
        self.R_core = np.hypot(corners[:, 0], corners[:, 1]).max()

    # ---- aggregate ----
    @property
    def all_pins(self):       return [p for a in self.assemblies for p in a.pins]
    @property
    def n_pins_total(self):   return self.N_assy * self.N_side ** 2
    @property
    def V_assemblies(self):   return self.N_assy * self.p_assy ** 2 * self.H
    @property
    def V_total(self):        return self.V_assemblies   # alias
    @property
    def V_core_cyl(self):     return np.pi * self.R_core ** 2 * self.H

    # ---- physics ----
    def sum_J0(self):
        return j0(c.J0_FIRST_ZERO * self._r / self.R_core).sum()

    def solve_H(self, P_target=c.P_THERMAL):
        """Closed-form H for target power. Mutates pins + self.H. Returns H."""
        A_fuel = self.assemblies[0].pins[0].A_fuel
        H = (P_target * np.pi) / (
            2 * c.N_FF * c.SIGMA_F * c.G_F_J * c.PHI_MAX * A_fuel * self.sum_J0()
        )
        self.H = H
        for a in self.assemblies:
            a.H = H
            for p in a.pins:
                p.H = H
        return H

    def total_power(self):
        """Sanity: sum per-pin power. Should match P_target after solve_H."""
        A_fuel = self.assemblies[0].pins[0].A_fuel
        return (c.N_FF * c.SIGMA_F * c.G_F_J * c.PHI_MAX
                * A_fuel * 2 * self.H / np.pi * self.sum_J0())

    def center_pin(self):
        idx = np.argmin(self._r)
        return self.all_pins[idx]

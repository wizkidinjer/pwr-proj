# core.py
import numpy as np
from assembly import Assembly


def _gauss_circle_layout(N_target):
    """
    Return N_target cell offsets (i, j) packed into a circle on a unit
    square grid, by including cells whose center lies within increasing
    radius R (Gauss circle problem). Raises if N_target is not achievable
    as an exact inclusion count.
    """
    for R in np.linspace(0, 10, 8000):
        rng = int(np.ceil(R)) + 1
        cells = [
            (i, j)
            for i in range(-rng, rng + 1)
            for j in range(-rng, rng + 1)
            if np.hypot(i, j) <= R
        ]
        if len(cells) == N_target:
            return sorted(cells)
    raise ValueError(f"No Gauss-circle layout of size {N_target}.")


# Assembly center offsets in units of assembly_pitch.
# Each pattern centered on origin. Small rubric-example patterns kept
# explicit; larger circle-filling patterns generated on demand.
LAYOUTS = {
    # -- rubric examples --
    5:  [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)],
    9:  [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)],
    21: (
        [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]
        + [(2, 0), (-2, 0), (0, 2), (0, -2)]
        + [(2, 1), (2, -1), (-2, 1), (-2, -1),
           (1, 2), (-1, 2), (1, -2), (-1, -2)]
    ),
    25: [(i, j) for i in range(-2, 3) for j in range(-2, 3)],

    # -- circle-filling (Gauss circle inclusion) --
    37: _gauss_circle_layout(37),   # fill 0.812
    57: _gauss_circle_layout(57),   # fill 0.806
    69: _gauss_circle_layout(69),   # fill 0.829
    89: _gauss_circle_layout(89),   # fill 0.872  best < 100
    97: _gauss_circle_layout(97),   # fill 0.846
}


class Core:
    def __init__(self, D_rod, assembly_pitch, N_assemblies,
                 N_side, clad_thick, H_active, MAX_P=15.5e6, P_DIFF=125e3):
        """
        D_rod:           rod outer diameter (m)
        assembly_pitch:  assembly edge length = center-to-center spacing (m)
        N_assemblies:    count of assemblies; must be a key in LAYOUTS
        N_side:          pins per edge inside each assembly
        clad_thick:      cladding thickness (m)
        H_active:        active fuel length (m)
        """
        if N_assemblies not in LAYOUTS:
            raise ValueError(
                f"N_assemblies={N_assemblies} not supported. "
                f"Pick from {sorted(LAYOUTS)}."
            )
        self.D_rod = D_rod
        self.assembly_pitch = assembly_pitch
        self.N_assemblies = N_assemblies
        self.N_side = N_side
        self.H = H_active
        self.MAX_P = MAX_P
        self.P_DIFF = P_DIFF

        self.assemblies = self._build_assemblies(clad_thick)
        self.R_core = self._compute_R_core()

    def _build_assemblies(self, clad_thick):
        offsets = LAYOUTS[self.N_assemblies]
        p = self.assembly_pitch
        asy = []
        for (i, j) in offsets:
            cx, cy = i * p, j * p
            asy.append(Assembly(
                N_side=self.N_side,
                assembly_pitch=p,
                center_xy=(cx, cy),
                pin_D=self.D_rod,
                clad_thick=clad_thick,
                H_active=self.H,
            ))
        return asy

    def _compute_R_core(self):
        """
        Circumscribed radius: smallest circle containing every assembly
        corner. That corner is the furthest fuel-region point from origin,
        so use it as R for flux shape J0(2.405 r / R).
        """
        r_max = 0.0
        for a in self.assemblies:
            xmin, ymin, xmax, ymax = a.bbox
            for (x, y) in [(xmin, ymin), (xmin, ymax),
                           (xmax, ymin), (xmax, ymax)]:
                r_max = max(r_max, np.hypot(x, y))
        return r_max

    @property
    def all_pins(self):
        """Flat list of every Pin in core."""
        return [p for a in self.assemblies for p in a.pins]

    @property
    def n_pins_total(self):
        return sum(a.n_pins for a in self.assemblies)

    @property
    def V_fuel_total(self):
        return sum(p.V_fuel for a in self.assemblies for p in a.pins)

    @property
    def V_clad_total(self):
        return sum(p.V_clad for a in self.assemblies for p in a.pins)

    @property
    def V_coolant_total(self):
        """Coolant volume inside assembly footprints only.
        Inter-assembly gaps = 0 (edge-tiled)."""
        return sum(a.V_coolant for a in self.assemblies)

    @property
    def V_core_cyl(self):
        """Volume of circumscribed cylinder (flux-domain reference)."""
        return np.pi * self.R_core**2 * self.H

    @property
    def V_assemblies(self):
        """
        Core volume = sum of assembly footprints x active height.
        This is the quantity to minimize per rubric.
            V = N_assemblies * assembly_pitch^2 * H
        """
        return self.N_assemblies * self.assembly_pitch**2 * self.H

    def solve_H_for_power(self, P_target, phi_max, N_ff, sigma_f, G_f):
        """
        Closed-form solve for active fuel length H that makes core deliver
        exactly P_target [W] under thin-rod flux assumption.

            H = P_target * pi
                / (2 * N_ff * sigma_f * G_f * phi_max * A_fuel * sum_J0)

        Mutates H on every pin and on self. Returns new H [m].
        """
        from scipy.special import j0
        A_fuel = self.assemblies[0].pins[0].A_fuel
        sum_J0 = sum(
            j0(2.405 * p.r / self.R_core)
            for a in self.assemblies for p in a.pins
        )
        H_new = (P_target * np.pi) / (
            2 * N_ff * sigma_f * G_f * phi_max * A_fuel * sum_J0
        )
        # mutate every pin + core
        for p in self.all_pins:
            p.H = H_new
        self.H = H_new
        return H_new

    def total_power(self, phi_max, N_ff, sigma_f, G_f):
        """Sanity-check: sum per-pin power. Should match P_target after solve."""
        return sum(
            p.power(self.R_core, phi_max, N_ff, sigma_f, G_f)
            for p in self.all_pins
        )

    def summary(self):
        return {
            "N_assemblies":    self.N_assemblies,
            "N_side":          self.N_side,
            "rods_per_assy":   self.N_side**2,
            "rods_total":      self.n_pins_total,
            "D_rod_mm":        self.D_rod * 1e3,
            "assy_pitch_mm":   self.assembly_pitch * 1e3,
            "rod_pitch_mm":    self.assemblies[0].pin_pitch * 1e3,
            "rod_spacing_mm":  self.assemblies[0].rod_spacing * 1e3,
            "P_over_D":        self.assemblies[0].pin_pitch / self.D_rod,
            "R_core_m":        self.R_core,
            "D_core_m":        2 * self.R_core,
            "H_core_m":        self.H,
            "V_core_m3":       self.V_core_cyl,
            "V_assy_m3":       self.V_assemblies,
            "V_fuel_m3":       self.V_fuel_total,
        }
# pin.py
import numpy as np
from scipy.special import j0


class Pin:
    def __init__(self, D_outer, clad_thick, H_active, loc_xy):
        """
        D_outer: fuel+clad outer diameter (m)
        clad_thick: cladding thickness (m)
        H_active: active fuel length (m)
        loc_xy: (x, y) pin center in core frame (m)
        """
        self.D_outer = D_outer
        self.clad_thick = clad_thick
        self.D_fuel = D_outer - 2 * clad_thick
        if self.D_fuel <= 0:
            raise ValueError(f"Clad too thick for D_outer={D_outer}")
        self.H = H_active
        self.x, self.y = loc_xy

    @property
    def r(self):
        return np.sqrt(self.x**2 + self.y**2)

    @property
    def A_fuel(self):
        return np.pi * (self.D_fuel / 2) ** 2

    @property
    def A_clad(self):
        return np.pi * ((self.D_outer / 2) ** 2 - (self.D_fuel / 2) ** 2)

    @property
    def V_fuel(self):
        return self.A_fuel * self.H

    @property
    def V_clad(self):
        return self.A_clad * self.H

    def flux_shape(self, R_ex, H_ex):
        """Dimensionless shape at pin, axially averaged over active H."""
        radial = j0(2.405 * self.r / R_ex)
        axial_avg = (2 * H_ex / (np.pi * self.H)) * np.sin(
            np.pi * self.H / (2 * H_ex)
        )
        return radial * axial_avg

    def radial_shape(self, R_core):
        """J0(2.405 r / R). Held constant across rod (thin-rod assumption)."""
        return j0(2.405 * self.r / R_core)

    def axial_integral(self):
        """
        Integral of cos(pi z / H) over z in [-H/2, H/2].
        Rubric: L_ex = H, so cos goes to zero at rod ends. Returns 2H/pi.
        """
        return 2 * self.H / np.pi

    def power(self, R_core, phi_max, N_ff, sigma_f, G_f):
        """
        Rod power [W] under thin-rod assumption:
            P = N_ff * sigma_f * G_f * phi_max
                * J0(2.405 r_pin / R_core) * A_fuel * (2 H / pi)
        """
        return (N_ff * sigma_f * G_f * phi_max
                * self.radial_shape(R_core)
                * self.A_fuel
                * self.axial_integral())
# assembly.py
from rod import Pin
import numpy as np
from config import MOD_RHO, MOD_MU, H_ACTIVE


class Assembly:
    def __init__(self, N_side, assembly_pitch, center_xy,
                 pin_D, clad_thick, H_active, MAX_P=None, P_DIFF=None):
        self.pin_D = pin_D
        self.N_side = N_side
        self.P_DIFF = P_DIFF
        self.MAX_P = MAX_P
        self.assembly_pitch = assembly_pitch
        self.P_outer = assembly_pitch
        self.pin_pitch = assembly_pitch / N_side
        if pin_D >= self.pin_pitch:
            raise ValueError(
                f"D={pin_D*1e3:.2f} mm >= rod pitch={self.pin_pitch*1e3:.2f} mm, pins overlap."
            )
        self.cx, self.cy = center_xy
        self.pins = self._build_pins(pin_D, clad_thick, H_active)

    def iterate_m_dot(self, initial_m_dot, verbose=False):
        mdot = initial_m_dot
        history = [mdot]
        for k in range(100):
            A_c = (self.pin_pitch**2 - np.pi * self.pin_D**2 / 4)
            um = mdot / 1000 / A_c
            HyD = 4 * A_c / (np.pi * self.pin_D)
            Re = MOD_RHO * um * HyD / MOD_MU
            a, b1, b2 = 0.1339, 0.09059, -0.09926
            C = a + b1 * (self.pin_pitch/self.pin_D - 1) + b2 * (self.pin_pitch/self.pin_D - 1)**2
            f = C / Re**0.18
            mdot_new = np.sqrt(self.MAX_P * 2 * MOD_RHO * A_c**2 / (f * H_ACTIVE / self.pin_D))
            print(mdot_new)
            history.append(mdot_new)
            if verbose:
                print(f"  iter {k:2d}: mdot={mdot_new:.4f}  Re={Re:.2e}  f={f:.4f}")
            if mdot_new - mdot < 0.01 * mdot:
                return mdot_new, history, True, k+1
            mdot = mdot_new
        return mdot, history, False, 100

    def _build_pins(self, D, clad, H):
        pins = []
        x0 = self.cx - self.P_outer / 2 + self.pin_pitch / 2
        y0 = self.cy - self.P_outer / 2 + self.pin_pitch / 2
        for i in range(self.N_side):
            for j in range(self.N_side):
                x = x0 + i * self.pin_pitch
                y = y0 + j * self.pin_pitch
                pins.append(Pin(D, clad, H, (x, y)))
        return pins

    @property
    def n_pins(self):
        return self.N_side ** 2

    @property
    def rod_spacing(self):
        return self.pin_pitch - self.pins[0].D_outer

    @property
    def V_coolant(self):
        H = self.pins[0].H
        V_tot = self.P_outer**2 * H
        return V_tot - sum(p.V_fuel + p.V_clad for p in self.pins)

    @property
    def bbox(self):
        half = self.P_outer / 2
        return (self.cx - half, self.cy - half, self.cx + half, self.cy + half)
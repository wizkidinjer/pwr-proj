# assembly.py
from pin import Pin


class Assembly:
    def __init__(self, N_side, pin_pitch, center_xy,
                 pin_D, clad_thick, H_active):
        """All lengths in m. N_side = pins per edge."""
        if pin_D >= pin_pitch:
            raise ValueError(f"D={pin_D} >= p={pin_pitch}, pins overlap")
        self.N_side = N_side
        self.pin_pitch = pin_pitch
        self.P_outer = N_side * pin_pitch
        self.cx, self.cy = center_xy
        self.pins = self._build_pins(pin_D, clad_thick, H_active)

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
    def V_coolant(self):
        H = self.pins[0].H
        V_tot = self.P_outer**2 * H
        return V_tot - sum(p.V_fuel + p.V_clad for p in self.pins)
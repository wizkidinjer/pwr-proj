# -*- coding: utf-8 -*-
from scipy.special import jv
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

'''REACTOR CONFIGURATION'''
EperF = 200 * 1.6022e-13  # joules per fission

class pwr:
    fuelrho = 8e20   # at/cm3
    sigf    = 585e-24  # cm2
    P       = 500e6
    f0      = 5e13

    def __init__(self, R, Z):
        self.R = R
        self.Z = Z

    def flux(self, r, z=0):
        """Flux at radial position r (cm) and axial position z (cm, default midplane)"""
        return self.f0 * jv(0, 2.405 * r / self.R) * np.cos(np.pi * z / self.Z)


# ── Reactor instance ──────────────────────────────────────────────────────────
R_ext = 211.8   # cm extrapolated radius
L_ext = 427.0   # cm extrapolated height
core  = pwr(R_ext, L_ext)

# ── 17x17 lattice geometry ────────────────────────────────────────────────────
N     = 17
pitch = 1.26    # cm rod pitch

offset = (N - 1) / 2
i_idx  = np.arange(N)
j_idx  = np.arange(N)
ii, jj = np.meshgrid(i_idx, j_idx)

x_cm = (ii - offset) * pitch   # physical x position of each rod (cm)
y_cm = (jj - offset) * pitch   # physical y position of each rod (cm)
r_cm = np.sqrt(x_cm**2 + y_cm**2)

# ── Flux at each rod (midplane z=0) ──────────────────────────────────────────
flux_grid = core.flux(r_cm, z=0)
# Rods outside extrapolated radius get zero (they don't exist)
flux_grid[r_cm >= R_ext] = 0.0

# ── 3D Bar Plot ───────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(11, 8))
ax  = fig.add_subplot(111, projection='3d')

# Flatten for bar3d
xpos   = ii.ravel().astype(float)
ypos   = jj.ravel().astype(float)
zpos   = np.zeros_like(xpos)
flux_f = flux_grid.ravel()

dx = dy = 0.8   # bar width (in lattice-index units)
dz = flux_f     # bar height = flux value

# Normalize flux for colormap
norm    = plt.Normalize(vmin=flux_f.min(), vmax=flux_f.max())
colors  = cm.plasma(norm(flux_f))

ax.bar3d(xpos - dx/2, ypos - dy/2, zpos,
         dx, dy, dz,
         color=colors, shade=True, zsort='average')

# ── Labels & formatting ───────────────────────────────────────────────────────
ax.set_xlabel('Rod column (i)', labelpad=8)
ax.set_ylabel('Rod row (j)', labelpad=8)
ax.set_zlabel('φ (n/cm²·s)', labelpad=8)
ax.set_title('17×17 PWR Lattice — Midplane Neutron Flux\n'
             r'$\phi(r) = \phi_0 \, J_0\!\left(\frac{2.405\,r}{R}\right)\cos\!\left(\frac{\pi z}{L}\right)$',
             pad=14)

ax.set_xticks(range(N))
ax.set_yticks(range(N))
ax.tick_params(axis='both', labelsize=6)

# Colorbar
mappable = cm.ScalarMappable(norm=norm, cmap='plasma')
mappable.set_array(flux_f)
cbar = fig.colorbar(mappable, ax=ax, shrink=0.5, pad=0.1, aspect=20)
cbar.set_label('Neutron flux (n/cm²·s)', fontsize=9)

ax.view_init(elev=30, azim=-60)
plt.tight_layout()
plt.savefig('pwr_lattice_flux.png', dpi=150)
plt.show()
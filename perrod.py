# -*- coding: utf-8 -*-
"""
v1_per_rod.py
Flux is defined over a single rod's local coordinates (r from rod axis, z along rod).
Every rod in the 17x17 lattice has the same flux distribution — same color gradient.
"""
from scipy.special import jv
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

EperF = 200 * 1.6022e-13

class pwr:
    fuelrho = 8e20
    sigf    = 585e-24
    P       = 500e6
    f0      = 5e13

    def __init__(self, R, Z):
        """
        R - rod radius (cm)
        Z - rod active length (cm)
        """
        self.R = R
        self.Z = Z

    def flux(self, r, z):
        """
        Flux over a single rod's local coordinates.
        r : radial distance from rod centerline (0 to R)
        z : axial position along rod (-Z/2 to Z/2)
        """
        return self.f0 * jv(0, 2.405 * r / self.R) * np.cos(np.pi * z / self.Z)


# ── Single rod geometry ───────────────────────────────────────────────────────
R_rod  = 0.4096   # cm  fuel pellet radius (Westinghouse)
L_rod  = 366.0    # cm  active fuel length

rod = pwr(R_rod, L_rod)

# Evaluate flux over the rod surface (r = R_rod, z sweeps full length)
n_z    = 40
n_th   = 24
z_vals = np.linspace(-L_rod / 2, L_rod / 2, n_z)
theta  = np.linspace(0, 2 * np.pi, n_th, endpoint=False)

# Surface flux at r = R_rod for each z
f_surface = rod.flux(R_rod, z_vals)          # shape (n_z,)
f_min, f_max = f_surface.min(), f_surface.max()
norm  = plt.Normalize(f_min, f_max)
cmap  = cm.plasma

# ── Build cylinder faces for one rod, colored by axial flux ──────────────────
def rod_faces_and_colors(xc, yc, z_offset=0):
    """
    Returns Poly3DCollection faces + facecolors for one rod.
    xc, yc    - rod center in the lattice (cm)
    z_offset  - base z position (all rods start at 0 here)
    """
    faces  = []
    colors = []
    xs = xc + R_rod * np.cos(theta)
    ys = yc + R_rod * np.sin(theta)

    for zi in range(n_z - 1):
        z0 = z_offset + z_vals[zi]
        z1 = z_offset + z_vals[zi+1]
        fval = (f_surface[zi] + f_surface[zi+1]) / 2
        col  = cmap(norm(fval))
        for k in range(n_th):
            k2 = (k + 1) % n_th
            face = [
                [xs[k],  ys[k],  z0],
                [xs[k2], ys[k2], z0],
                [xs[k2], ys[k2], z1],
                [xs[k],  ys[k],  z1],
            ]
            faces.append(face)
            colors.append(col)

    return faces, colors


# ── 17x17 lattice ─────────────────────────────────────────────────────────────
N      = 17
pitch  = 1.26   # cm
offset = (N - 1) / 2

fig = plt.figure(figsize=(11, 9))
ax  = fig.add_subplot(111, projection='3d')

all_faces  = []
all_colors = []

for i in range(N):
    for j in range(N):
        xc = (i - offset) * pitch
        yc = (j - offset) * pitch
        faces, colors = rod_faces_and_colors(xc, yc)
        all_faces.extend(faces)
        all_colors.extend(colors)

poly = Poly3DCollection(all_faces, facecolors=all_colors, edgecolor='none', alpha=0.95)
ax.add_collection3d(poly)

# ── Axes & labels ─────────────────────────────────────────────────────────────
half = (N / 2) * pitch + pitch
ax.set_xlim(-half, half)
ax.set_ylim(-half, half)
ax.set_zlim(-L_rod / 2, L_rod / 2)

ax.set_xlabel('x (cm)', labelpad=6)
ax.set_ylabel('y (cm)', labelpad=6)
ax.set_zlabel('z (cm)', labelpad=6)
ax.set_title('17×17 PWR Lattice — Per-rod flux distribution\n'
             r'$\phi(r,z)=\phi_0\,J_0\!\left(\frac{2.405\,r}{R_{rod}}\right)'
             r'\cos\!\left(\frac{\pi z}{L}\right)$  (same on every rod)',
             pad=10)

mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
mappable.set_array(f_surface)
cbar = fig.colorbar(mappable, ax=ax, shrink=0.45, pad=0.08)
cbar.set_label('φ (n/cm²·s)')

ax.view_init(elev=25, azim=-55)
plt.tight_layout()
plt.savefig('v1_per_rod.png', dpi=150)
plt.show()
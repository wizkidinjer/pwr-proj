# -*- coding: utf-8 -*-
"""
v2_continuous_field.py
Flux defined as a continuous field phi(r, z) over the whole core.
Rods are cylinders colored by sampling that field at each rod's (r, z).
Each rod's color gradient differs — rods near the center are brighter overall
because the Bessel radial term scales them, AND all share the cosine axial shape.
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
        R - extrapolated core radius (cm)
        Z - extrapolated core height (cm)
        """
        self.R = R
        self.Z = Z

    def flux(self, r, z):
        """
        Continuous flux field over the core.
        r : radial distance from core centerline (cm)
        z : axial position from core midplane (cm), z=0 at center
        """
        return self.f0 * jv(0, 2.405 * r / self.R) * np.cos(np.pi * z / self.Z)


# ── Core instance ─────────────────────────────────────────────────────────────
R_ext = 211.8   # cm extrapolated radius
L_ext = 427.0   # cm extrapolated height
core  = pwr(R_ext, L_ext)

# ── Rod geometry ──────────────────────────────────────────────────────────────
N      = 17
pitch  = 1.26     # cm
offset = (N - 1) / 2
R_rod  = 0.4096   # cm fuel pellet radius
L_rod  = 366.0    # cm active fuel length

n_z  = 40
n_th = 24
z_vals = np.linspace(-L_rod / 2, L_rod / 2, n_z)
theta  = np.linspace(0, 2 * np.pi, n_th, endpoint=False)

# ── Global flux range (for shared colormap) ───────────────────────────────────
# Sample flux at core center rod (r~0) and edge rod across all z
r_center = 0.0
r_edge   = offset * pitch  # outermost rod radius from core center
f_peak   = core.flux(r_center, 0.0)
f_low    = core.flux(r_edge,   L_rod / 2)
norm     = plt.Normalize(f_low, f_peak)
cmap     = cm.plasma

# ── Build cylinder faces colored by field sample at (r_rod, z) ───────────────
def rod_faces_and_colors(xc, yc):
    r_rod = np.sqrt(xc**2 + yc**2)   # rod's radial position in core (cm)
    faces  = []
    colors = []
    xs = xc + R_rod * np.cos(theta)
    ys = yc + R_rod * np.sin(theta)

    for zi in range(n_z - 1):
        z_mid = (z_vals[zi] + z_vals[zi + 1]) / 2
        fval  = core.flux(r_rod, z_mid)   # <-- field sampled at this rod's (r, z)
        col   = cmap(norm(fval))
        z0, z1 = z_vals[zi], z_vals[zi + 1]
        for k in range(n_th):
            k2 = (k + 1) % n_th
            faces.append([
                [xs[k],  ys[k],  z0],
                [xs[k2], ys[k2], z0],
                [xs[k2], ys[k2], z1],
                [xs[k],  ys[k],  z1],
            ])
            colors.append(col)

    return faces, colors

# ── Plot ──────────────────────────────────────────────────────────────────────
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

half = (N / 2) * pitch + pitch
ax.set_xlim(-half, half)
ax.set_ylim(-half, half)
ax.set_zlim(-L_rod / 2, L_rod / 2)

ax.set_xlabel('x (cm)', labelpad=6)
ax.set_ylabel('y (cm)', labelpad=6)
ax.set_zlabel('z (cm)', labelpad=6)
ax.set_title('17×17 PWR Lattice — Continuous core flux field\n'
             r'$\phi(r,z)=\phi_0\,J_0\!\left(\frac{2.405\,r}{R_{core}}\right)'
             r'\cos\!\left(\frac{\pi z}{L}\right)$  sampled at each rod $(r,z)$',
             pad=10)

mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
mappable.set_array([f_low, f_peak])
cbar = fig.colorbar(mappable, ax=ax, shrink=0.45, pad=0.08)
cbar.set_label('φ (n/cm²·s)')

ax.view_init(elev=25, azim=-55)
plt.tight_layout()
plt.savefig('v2_continuous_field.png', dpi=150)
plt.show()
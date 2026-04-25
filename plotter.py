# plotter.py
"""Top-down x-y core cross section. Optional, not required by main flow."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
from scipy.special import j0
import config as c


def plot_core_xy(core, flux_contour=True, show_rods=True, savepath=None):
    fig, ax = plt.subplots(figsize=(9, 9))
    R = core.R_core
    if flux_contour:
        N = 300
        xs = np.linspace(-R, R, N); X, Y = np.meshgrid(xs, xs)
        rr = np.hypot(X, Y)
        phi = np.where(rr <= R, j0(c.J0_FIRST_ZERO * rr / R), np.nan)
        cf = ax.contourf(X, Y, phi, levels=15, cmap="viridis", alpha=0.55)
        plt.colorbar(cf, ax=ax, label=r"$\phi/\phi_{max}$")
    for a in core.assemblies:
        x0, y0, x1, y1 = a.bbox
        ax.add_patch(Rectangle((x0, y0), x1-x0, y1-y0, fill=False,
                               edgecolor="white", linewidth=1.4, zorder=2))
        ax.plot(a.cx, a.cy, "w+", ms=6, zorder=3)
    if show_rods:
        rad = core.D_rod / 2
        for p in core.all_pins:
            ax.add_patch(Circle((p.x, p.y), rad, facecolor="black",
                                edgecolor="none", zorder=4))
    ax.add_patch(Circle((0, 0), R, fill=False, edgecolor="red",
                        linewidth=2.2, linestyle="--", zorder=5,
                        label=f"R_core={R*1e3:.1f} mm"))
    pad = 0.08*R
    ax.set_xlim(-R-pad, R+pad); ax.set_ylim(-R-pad, R+pad)
    ax.set_aspect("equal"); ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)")
    ax.set_title(f"{core.N_assy} assys | {core.n_pins_total} rods | D_core={2*R:.3f} m")
    ax.legend(loc="upper right"); ax.grid(alpha=0.25)
    if savepath: fig.savefig(savepath, dpi=150, bbox_inches="tight")
    return fig, ax


def plot_core_3d(core, show_rods=True, savepath=None):
    """3D wireframe of core: assembly boxes + flux-cylinder bounds.
       show_rods: render every rod as a thin cylinder (slow for big cores)."""
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    fig = plt.figure(figsize=(10, 9))
    ax  = fig.add_subplot(111, projection="3d")
    H, R = core.H, core.R_core
    z0, z1 = -H/2, H/2

    # ---- assembly bounding boxes (white edges, transparent faces) ----
    for a in core.assemblies:
        x0, y0, x1, y1 = a.bbox
        verts = np.array([
            [[x0,y0,z0],[x1,y0,z0],[x1,y1,z0],[x0,y1,z0]],
            [[x0,y0,z1],[x1,y0,z1],[x1,y1,z1],[x0,y1,z1]],
            [[x0,y0,z0],[x1,y0,z0],[x1,y0,z1],[x0,y0,z1]],
            [[x1,y0,z0],[x1,y1,z0],[x1,y1,z1],[x1,y0,z1]],
            [[x1,y1,z0],[x0,y1,z0],[x0,y1,z1],[x1,y1,z1]],
            [[x0,y1,z0],[x0,y0,z0],[x0,y0,z1],[x0,y1,z1]],
        ])
        ax.add_collection3d(Poly3DCollection(
            verts, facecolors=(0.4,0.6,0.9,0.08),
            edgecolors="steelblue", linewidths=0.8))

    # ---- flux-cylinder wireframe (red dashed) ----
    th = np.linspace(0, 2*np.pi, 80)
    cx, cy = R*np.cos(th), R*np.sin(th)
    ax.plot(cx, cy, z0, "r--", lw=1.5)
    ax.plot(cx, cy, z1, "r--", lw=1.5)
    for k in range(0, 80, 10):
        ax.plot([cx[k]]*2, [cy[k]]*2, [z0, z1], "r--", lw=0.8, alpha=0.6)

    # ---- rods (optional) ----
    if show_rods:
        rad = core.D_rod / 2
        th_r = np.linspace(0, 2*np.pi, 12)
        for p in core.all_pins:
            xs = p.x + rad*np.cos(th_r); ys = p.y + rad*np.sin(th_r)
            for zz in (z0, z1):
                ax.plot(xs, ys, zz, "k-", lw=0.4, alpha=0.7)
            ax.plot([p.x, p.x], [p.y, p.y], [z0, z1], "k-", lw=0.3, alpha=0.5)

    ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)"); ax.set_zlabel("z (m)")
    ax.set_title(f"3D core: {core.N_assy} assys | {core.n_pins_total} rods | "
                 f"H={H:.2f} m | D_core={2*R:.2f} m")
    pad = 0.05*R
    ax.set_xlim(-R-pad, R+pad); ax.set_ylim(-R-pad, R+pad); ax.set_zlim(z0, z1)
    ax.set_box_aspect((2*(R+pad), 2*(R+pad), H))   # true x:y:z scale
    if savepath: fig.savefig(savepath, dpi=150, bbox_inches="tight")
    return fig, ax
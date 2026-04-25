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

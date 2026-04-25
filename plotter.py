# plotter.py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle


def plot_core_xy(core, flux_contour=True, show_rods=True,
                 rod_stride=1, savepath=None):
    """
    Top-down (x-y) axial cross-section of core.

    core:         Core instance
    flux_contour: overlay J0(2.405 r / R_core) contours
    show_rods:    draw individual rod circles (slow for 7k+ rods)
    rod_stride:   draw every Nth rod if show_rods (1 = all)
    savepath:     if given, save to path
    """
    fig, ax = plt.subplots(figsize=(9, 9))

    # --- flux contour backdrop ---
    if flux_contour:
        R = core.R_core
        N = 300
        xs = np.linspace(-R, R, N)
        ys = np.linspace(-R, R, N)
        X, Y = np.meshgrid(xs, ys)
        r = np.hypot(X, Y)
        from scipy.special import j0
        phi = np.where(r <= R, j0(2.405 * r / R), np.nan)
        cf = ax.contourf(X, Y, phi, levels=15, cmap="viridis", alpha=0.55)
        plt.colorbar(cf, ax=ax, label=r"$\phi / \phi_{max}$  (radial, J$_0$)")

    # --- assembly bounding boxes ---
    for a in core.assemblies:
        xmin, ymin, xmax, ymax = a.bbox
        ax.add_patch(Rectangle(
            (xmin, ymin), xmax - xmin, ymax - ymin,
            fill=False, edgecolor="white", linewidth=1.4, zorder=2,
        ))
        ax.plot(a.cx, a.cy, "w+", markersize=6, zorder=3)

    # --- rods ---
    if show_rods:
        r_rod = core.D_rod / 2
        for k, p in enumerate(core.all_pins):
            if k % rod_stride != 0:
                continue
            ax.add_patch(Circle(
                (p.x, p.y), r_rod,
                facecolor="black", edgecolor="none", zorder=4,
            ))

    # --- core circumscribed circle (flux outer radius R) ---
    ax.add_patch(Circle(
        (0, 0), core.R_core,
        fill=False, edgecolor="red", linewidth=2.2,
        linestyle="--", zorder=5, label=f"R_core = {core.R_core*1e3:.1f} mm",
    ))

    # --- cosmetics ---
    pad = 0.08 * core.R_core
    ax.set_xlim(-core.R_core - pad, core.R_core + pad)
    ax.set_ylim(-core.R_core - pad, core.R_core + pad)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(
        f"Core cross-section  |  {core.N_assemblies} assemblies  "
        f"|  {core.n_pins_total} rods  |  D_core = {2*core.R_core:.3f} m"
    )
    ax.legend(loc="upper right")
    ax.grid(alpha=0.25)

    if savepath:
        fig.savefig(savepath, dpi=150, bbox_inches="tight")
    return fig, ax


if __name__ == "__main__":
    # quick demo: pull from config, build core, plot
    import config as c
    from core import Core

    core = Core(
        D_rod=c.D_ROD_OUTER,
        assembly_pitch=c.N_SIDE * c.PIN_PITCH,   # same as old P_outer
        N_assemblies=c.N_ASSEMBLIES,
        N_side=c.N_SIDE,
        clad_thick=c.CLAD_THICK,
        H_active=c.H_ACTIVE,
    )
    print("Summary:")
    for k, v in core.summary().items():
        print(f"  {k:16s} = {v}")

    plot_core_xy(core, savepath="core_xy.png")
    plt.show()
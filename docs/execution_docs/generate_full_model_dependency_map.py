from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)


PALETTE = {
    "ink": "#0f172a",
    "muted": "#475569",
    "lane_exec": "#eef2ff",
    "lane_shared": "#ecfeff",
    "lane_branch": "#f8fafc",
    "lane_output": "#f8fafc",
    "input": "#e8eefc",
    "exec": "#dbeafe",
    "shared": "#d1fae5",
    "shared_alt": "#fef3c7",
    "heat": "#fee2e2",
    "cool": "#dcfce7",
    "annual": "#fde68a",
    "output": "#e5e7eb",
}


def fig_ax():
    fig, ax = plt.subplots(figsize=(28, 20))
    ax.set_xlim(0, 220)
    ax.set_ylim(0, 165)
    ax.axis("off")
    fig.suptitle(
        "07. Current Full-Model Dependency Map",
        fontsize=24,
        fontweight="bold",
        y=0.988,
    )
    ax.text(
        110,
        158.5,
        "Python entrypoint + conditional MATLAB refresh/export + shared physics kernels + heating, cooling, and annual branches",
        ha="center",
        va="center",
        fontsize=12.5,
        color=PALETTE["muted"],
    )
    return fig, ax


def rounded_box(ax, x, y, w, h, fc, ec="#1f2937", lw=1.9, radius=1.8, alpha=1.0):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.45,rounding_size={radius}",
        linewidth=lw,
        edgecolor=ec,
        facecolor=fc,
        alpha=alpha,
        zorder=1.0,
    )
    ax.add_patch(patch)
    return patch


def lane(ax, x, y, w, h, title, fc):
    rounded_box(ax, x, y, w, h, fc=fc, ec="#cbd5e1", lw=1.4, radius=2.3)
    if title:
        ax.text(x + 2.5, y + h - 2.0, title, ha="left", va="top", fontsize=14, fontweight="bold", color=PALETTE["ink"])


def panel(ax, x, y, w, h, title, fc, ec):
    rounded_box(ax, x, y, w, h, fc=fc, ec=ec, lw=1.7, radius=2.0, alpha=0.38)
    rounded_box(ax, x + 2.0, y + h - 5.4, w - 4.0, 4.3, fc="#ffffff", ec=ec, lw=1.3, radius=1.2, alpha=0.98)
    ax.text(x + w / 2, y + h - 3.2, title, ha="center", va="center", fontsize=9.8, fontweight="bold", color=ec)


def box(ax, x, y, w, h, title, body, fc, ec="#1f2937", title_fs=11.5, body_fs=9.8):
    rounded_box(ax, x, y, w, h, fc=fc, ec=ec, lw=1.8, radius=1.6)
    ax.text(x + w / 2, y + h - 1.55, title, ha="center", va="top", fontsize=title_fs, fontweight="bold", color=PALETTE["ink"])
    ax.text(
        x + w / 2,
        y + h * 0.38,
        body,
        ha="center",
        va="center",
        fontsize=body_fs,
        color=PALETTE["ink"],
        linespacing=1.15,
    )


def arrow(ax, p0, p1, text=None, color="#111827", lw=2.0, ms=16, rad=0.0, text_xy=None, fs=9.5):
    patch = FancyArrowPatch(
        p0,
        p1,
        arrowstyle="-|>",
        mutation_scale=ms,
        linewidth=lw,
        color=color,
        connectionstyle=f"arc3,rad={rad}",
        zorder=2.0,
    )
    ax.add_patch(patch)
    if text:
        tx, ty = text_xy if text_xy is not None else ((p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2)
        ax.text(tx, ty, text, ha="center", va="center", fontsize=fs, color=color, zorder=4)


def note(ax, x, y, text, fs=10.0, color=None, ha="left"):
    ax.text(x, y, text, ha=ha, va="center", fontsize=fs, color=color or PALETTE["muted"])


def build_map():
    fig, ax = fig_ax()

    lane(ax, 2, 126, 216, 28, "Execution, configuration, and external data", PALETTE["lane_exec"])
    lane(ax, 2, 86, 216, 34, "Shared physics stack used across the full workflow", PALETTE["lane_shared"])
    lane(ax, 2, 16, 216, 64, "", PALETTE["lane_branch"])
    lane(ax, 2, 2, 216, 12, "", PALETTE["lane_output"])

    # Top lane: execution / inputs
    box(
        ax,
        4,
        131,
        28,
        18,
        "Study config + CLI",
        "run modes\nmodel overrides\nparallel settings\nactive vessel config",
        fc=PALETTE["input"],
    )
    box(
        ax,
        38,
        131,
        28,
        18,
        "heat_transfer_study.py",
        "load config\nsync aliases\nresolve paths\ncreate output tree",
        fc=PALETTE["exec"],
    )
    box(
        ax,
        72,
        131,
        28,
        18,
        "Mode dispatch",
        "heating branch\ncooling branch\nyear-round branch\nroot summary index",
        fc=PALETTE["exec"],
    )
    box(
        ax,
        106,
        131,
        28,
        18,
        "Conditional refresh/export",
        "compare override signature\nreuse compatible export\nor trigger MATLAB refresh",
        fc=PALETTE["exec"],
        body_fs=9.3,
    )
    box(
        ax,
        140,
        131,
        28,
        18,
        "analyze_vermicomposter_\nheater.py",
        "build temp export script\nlaunch matlab -batch\nannotate payload JSON",
        fc=PALETTE["exec"],
        body_fs=9.3,
    )
    box(
        ax,
        174,
        131,
        40,
        18,
        "External sources",
        "NOAA climate + events\nEIA electricity price\nspot-cooler and blower specs",
        fc="#e0f2fe",
        body_fs=9.5,
    )

    # Shared kernels
    box(
        ax,
        6,
        92,
        36,
        22,
        "Geometry + topology",
        "vessel configuration\nwire length + electrical strings\nflow split\nwetted area + tube layout",
        fc=PALETTE["shared_alt"],
        body_fs=9.4,
    )
    box(
        ax,
        48,
        92,
        36,
        22,
        "Air + moisture properties",
        "ideal gas + Sutherland\nhumidity ratio + dew point\nIAPWS saturation\nlatent heat",
        fc=PALETTE["shared"],
        body_fs=9.4,
    )
    box(
        ax,
        90,
        92,
        36,
        22,
        "Correlations + losses",
        "Churchill-Bernstein\nHausen / Gnielinski\nChurchill friction\nIdelchik losses",
        fc=PALETTE["shared"],
        body_fs=9.2,
    )
    box(
        ax,
        132,
        92,
        36,
        22,
        "Core reduced-order solvers",
        "heater tube\nperforated aeration tube\nevaporation limits\ntwo-node bed + wall UA",
        fc=PALETTE["shared"],
        body_fs=9.2,
    )
    box(
        ax,
        174,
        92,
        36,
        22,
        "Post-processing kernels",
        "radial plume cooling\nhabitat exclusion\ncross-section layout\ntrade-study plots",
        fc=PALETTE["shared_alt"],
        body_fs=9.3,
    )
    # Branch panels
    panel(ax, 4, 20, 68, 56, "Heating branch\nMATLAB candidates + Python recommendation", fc=PALETTE["heat"], ec="#b91c1c")
    panel(ax, 76, 20, 68, 56, "Cooling branch\nPython summer workflow", fc=PALETTE["cool"], ec="#047857")
    panel(ax, 148, 20, 68, 56, "Annual branch\nclimate baseline + daily re-solve", fc=PALETTE["annual"], ec="#a16207")

    # Heating branch boxes
    box(
        ax,
        10,
        58.0,
        56,
        9.6,
        "MATLAB operating sweep",
        "voltage x flow x geometry\nfor each active vessel configuration",
        fc="#fef2f2",
        title_fs=10.4,
        body_fs=8.1,
    )
    box(
        ax,
        10,
        46.2,
        56,
        9.6,
        "Per-point heating physics",
        "heater tube + aeration tube\npressure drop, water loss,\ntwo-node bed temperatures",
        fc="#fff1f2",
        title_fs=10.4,
        body_fs=8.0,
    )
    box(
        ax,
        10,
        34.4,
        56,
        9.6,
        "Discrete candidate table",
        "export JSON\nconfiguration comparison\nwire / geometry trade sweeps",
        fc="#fff1f2",
        title_fs=10.4,
        body_fs=8.0,
    )
    box(
        ax,
        10,
        22.6,
        56,
        9.6,
        "Python heating selection",
        "effective limits + metrics\nhard constraints\nNSGA-II Pareto + report priority",
        fc="#ffe4e6",
        title_fs=10.4,
        body_fs=8.0,
    )

    # Cooling branch boxes
    box(
        ax,
        82,
        58.0,
        56,
        9.6,
        "Spot cooler + assist blower",
        "supply state\ndehumidification / condensate\navailable static pressure",
        fc="#ecfdf5",
        title_fs=10.4,
        body_fs=8.0,
    )
    box(
        ax,
        82,
        46.2,
        56,
        9.6,
        "Summer cooling sweep",
        "flow array + bed-reference bisection\nPython aeration / evaporation / dP",
        fc="#ecfdf5",
        title_fs=10.4,
        body_fs=8.0,
    )
    box(
        ax,
        82,
        34.4,
        56,
        9.6,
        "Cooling metrics",
        "summer hard-safe residuals\nnegative Q_to_bed required\ncooling capacity + shortfall",
        fc="#ecfdf5",
        title_fs=10.4,
        body_fs=8.0,
    )
    box(
        ax,
        82,
        22.6,
        56,
        9.6,
        "Python cooling selection",
        "NSGA-II Pareto\nreport priority\nbest-available fallback",
        fc="#d1fae5",
        title_fs=10.4,
        body_fs=8.1,
    )

    # Annual branch boxes
    box(
        ax,
        154,
        58.0,
        56,
        9.6,
        "Representative climate year",
        "periodic Gaussian kernel baseline\nfreeze and heat-wave overrides",
        fc="#fffbeb",
        title_fs=10.4,
        body_fs=8.0,
    )
    box(
        ax,
        154,
        46.2,
        56,
        9.6,
        "Selected design points",
        "take the reported heating and cooling\nreference points from the branch payloads",
        fc="#fffbeb",
        title_fs=10.4,
        body_fs=7.9,
    )
    box(
        ax,
        154,
        34.4,
        56,
        9.6,
        "Daily ambient-specific re-solve",
        "recompute capacity, power,\nwater, and duty at each day",
        fc="#fffbeb",
        title_fs=10.4,
        body_fs=8.0,
    )
    box(
        ax,
        154,
        22.6,
        56,
        9.6,
        "Annual totals",
        "electricity\noperating cost\nwater replacement\nunmet-load equivalent",
        fc="#fef3c7",
        title_fs=10.4,
        body_fs=8.0,
    )

    # Output boxes
    box(
        ax,
        10,
        3.6,
        56,
        8.8,
        "Heating outputs",
        "constraint maps, trade studies,\nPareto front, summary text",
        fc=PALETTE["output"],
        title_fs=10.4,
        body_fs=7.7,
    )
    box(
        ax,
        82,
        3.6,
        56,
        8.8,
        "Cooling outputs",
        "spot-cooler curves, Pareto front,\nradial profile, summary text",
        fc=PALETTE["output"],
        title_fs=10.4,
        body_fs=7.7,
    )
    box(
        ax,
        154,
        3.6,
        56,
        8.8,
        "Year-round outputs",
        "energy / cost plot\nannual summary text",
        fc=PALETTE["output"],
        title_fs=10.4,
        body_fs=7.7,
    )
    note(ax, 110, 1.1, "Heating recommendation is Python-side, cooling is fully Python-side, and annual totals consume the reported branch design points.", fs=8.8, ha="center")

    # Main execution arrows
    arrow(ax, (32, 140), (38, 140), color="#1d4ed8")
    arrow(ax, (66, 140), (72, 140), color="#1d4ed8")
    arrow(ax, (100, 140), (106, 140), color="#1d4ed8")
    arrow(ax, (134, 140), (140, 140), color="#1d4ed8")

    # Mode dispatch arrows down
    arrow(ax, (86, 131), (38, 69.5), color="#2563eb", rad=0.08)
    arrow(ax, (86, 131), (110, 69.5), color="#16a34a", rad=-0.03)
    arrow(ax, (86, 131), (182, 69.5), color="#a16207", rad=-0.10)

    # Refresh/export path to heating sweep
    arrow(ax, (154, 131), (38, 69.0), color="#b91c1c", rad=0.12)

    # External data arrows
    arrow(ax, (194, 131), (110, 69.0), color="#0f766e", rad=0.08)
    arrow(ax, (194, 131), (182, 69.0), color="#a16207", rad=0.02)

    # Shared kernel flow
    arrow(ax, (42, 103), (48, 103), color="#64748b")
    arrow(ax, (84, 103), (90, 103), color="#64748b")
    arrow(ax, (126, 103), (132, 103), color="#64748b")
    arrow(ax, (168, 103), (174, 103), color="#64748b")
    arrow(ax, (24, 92), (26, 69.0), color="#b45309", rad=-0.04, lw=1.8, ms=14)
    arrow(ax, (66, 92), (38, 69.0), color="#0f766e", rad=0.16, lw=1.8, ms=14)
    arrow(ax, (108, 92), (54, 69.0), color="#0f766e", rad=0.10, lw=1.8, ms=14)
    arrow(ax, (150, 92), (66, 69.0), color="#0f766e", rad=0.12, lw=1.8, ms=14)
    arrow(ax, (24, 92), (96, 69.0), color="#b45309", rad=0.08, lw=1.8, ms=14)
    arrow(ax, (66, 92), (110, 69.0), color="#0f766e", rad=-0.04, lw=1.8, ms=14)
    arrow(ax, (108, 92), (122, 69.0), color="#0f766e", lw=1.8, ms=14)
    arrow(ax, (150, 92), (138, 69.0), color="#0f766e", rad=0.04, lw=1.8, ms=14)
    arrow(ax, (150, 92), (182, 69.0), color="#a16207", rad=-0.14)
    arrow(ax, (192, 92), (38, 12.6), color="#15803d", rad=0.26, lw=1.8, ms=14)
    arrow(ax, (192, 92), (110, 12.6), color="#15803d", rad=0.18, lw=1.8, ms=14)

    # Vertical branch arrows
    for x in (38, 110, 182):
        arrow(ax, (x, 58.0), (x, 55.8), color=PALETTE["ink"], ms=14)
        arrow(ax, (x, 46.2), (x, 44.0), color=PALETTE["ink"], ms=14)
        arrow(ax, (x, 34.4), (x, 32.2), color=PALETTE["ink"], ms=14)

    # Cross-branch annual dependencies
    arrow(ax, (66, 27.4), (182, 69.0), color="#a16207", rad=-0.38)
    arrow(ax, (138, 27.4), (182, 69.0), color="#a16207", rad=-0.20)

    # Outputs
    arrow(ax, (38, 22.6), (38, 12.4), color="#b91c1c", rad=0.0)
    arrow(ax, (110, 22.6), (110, 12.4), color="#047857", rad=0.0)
    arrow(ax, (182, 22.6), (182, 12.4), color="#a16207", rad=0.0)

    svg_path = OUT_DIR / "full_model_dependency_map_current.svg"
    png_path = OUT_DIR / "full_model_dependency_map_current.png"
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=220, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved {svg_path}")
    print(f"Saved {png_path}")


if __name__ == "__main__":
    build_map()

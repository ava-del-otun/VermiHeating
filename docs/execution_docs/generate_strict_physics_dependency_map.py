from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)


PALETTE = {
    "ink": "#0f172a",
    "muted": "#475569",
    "lane_input": "#eef2ff",
    "lane_shared": "#ecfeff",
    "lane_heat": "#fff1f2",
    "lane_cool": "#ecfdf5",
    "lane_annual": "#fffbeb",
    "input_ambient": "#dbeafe",
    "input_geom": "#fef3c7",
    "input_branch": "#e9d5ff",
    "input_heat": "#fee2e2",
    "input_cool": "#d1fae5",
    "shared_a": "#dbeafe",
    "shared_b": "#dcfce7",
    "shared_c": "#fef3c7",
    "shared_d": "#fae8ff",
    "shared_e": "#fee2e2",
    "heat": "#ffe4e6",
    "cool": "#dcfce7",
    "annual": "#fde68a",
    "post": "#e5e7eb",
    "blue": "#2563eb",
    "teal": "#0f766e",
    "green": "#15803d",
    "orange": "#b45309",
    "red": "#b91c1c",
    "gold": "#a16207",
    "gray": "#64748b",
}


def fig_ax():
    fig, ax = plt.subplots(figsize=(34, 26))
    ax.set_xlim(0, 270)
    ax.set_ylim(0, 246)
    ax.axis("off")
    fig.suptitle(
        "08. Strict Physics Dependency Map",
        fontsize=24,
        fontweight="bold",
        y=0.988,
    )
    ax.text(
        135,
        238.5,
        "Grouped variable/equation nodes for the live heating, cooling, annual, and post-processing physics",
        ha="center",
        va="center",
        fontsize=12.5,
        color=PALETTE["muted"],
    )
    return fig, ax


def rounded_box(ax, x, y, w, h, fc, ec="#1f2937", lw=1.8, radius=1.8, alpha=1.0, z=1.0):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.45,rounding_size={radius}",
        linewidth=lw,
        edgecolor=ec,
        facecolor=fc,
        alpha=alpha,
        zorder=z,
    )
    ax.add_patch(patch)
    return patch


def lane(ax, x, y, w, h, title, fc):
    rounded_box(ax, x, y, w, h, fc=fc, ec="#cbd5e1", lw=1.4, radius=2.3, z=0.3)
    ax.text(x + 3, y + h - 2.2, title, ha="left", va="top", fontsize=13.2, fontweight="bold", color=PALETTE["ink"])


def box(ax, x, y, w, h, title, body, fc, ec="#1f2937", title_fs=10.5, body_fs=7.7):
    rounded_box(ax, x, y, w, h, fc=fc, ec=ec, lw=1.5, radius=1.5, z=1.0)
    ax.text(x + w / 2, y + h - 1.45, title, ha="center", va="top", fontsize=title_fs, fontweight="bold", color=PALETTE["ink"], zorder=3)
    ax.text(
        x + w / 2,
        y + h * 0.40,
        body,
        ha="center",
        va="center",
        fontsize=body_fs,
        color=PALETTE["ink"],
        linespacing=1.12,
        zorder=3,
    )


def arrow(ax, p0, p1, color="#111827", lw=1.8, ms=14, rad=0.0):
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


def center_top(x, y, w, h):
    return (x + w / 2, y + h)


def center_bottom(x, y, w, h):
    return (x + w / 2, y)


def left_mid(x, y, w, h):
    return (x, y + h / 2)


def right_mid(x, y, w, h):
    return (x + w, y + h / 2)


def build_map():
    fig, ax = fig_ax()

    lane(ax, 2, 194, 266, 38, "Primary inputs and prescribed drivers", PALETTE["lane_input"])
    lane(ax, 2, 150, 266, 34, "Shared derived variables and governing closures", PALETTE["lane_shared"])
    lane(ax, 2, 102, 266, 38, "Heating physics branch", PALETTE["lane_heat"])
    lane(ax, 2, 56, 266, 38, "Cooling physics branch", PALETTE["lane_cool"])
    lane(ax, 2, 12, 266, 38, "Annual re-solve and post-processing layer", PALETTE["lane_annual"])

    x_cols = [4, 58, 112, 166, 220]
    w = 46

    # Inputs
    box(ax, x_cols[0], 200, w, 24, "Ambient + climate\ndrivers", "T_inf, p, RH_in\nmonthly mean/std\nfreeze + heat-wave stats", PALETTE["input_ambient"], body_fs=7.5)
    box(ax, x_cols[1], 200, w, 24, "Bin + wall\ninputs", "L, W, H, fill fraction\nk_bin, rho_bulk, c_p\nsheet/insulation\ntop + vent geometry", PALETTE["input_geom"], body_fs=7.2)
    box(ax, x_cols[2], 200, w, 24, "Aeration +\nperforation inputs", "N_tubes, L_aer, ID/OD\nrelease fraction\nholes/tube, d_h\nheader/splitter/connector", PALETTE["input_branch"], body_fs=7.0)
    box(ax, x_cols[3], 200, w, 24, "Heater + wire +\nsupply inputs", "V, string counts\nD_wire, D_coil, pitch/span\nrho_wire, alpha", PALETTE["input_heat"], body_fs=7.5)
    box(ax, x_cols[4], 200, w, 24, "Cooling equipment\ninputs", "T_supply,target, COP\nQ_max,cool\nrated condensate\nblower rated/shutoff curve", PALETTE["input_cool"], body_fs=7.0)

    # Shared closures
    box(ax, x_cols[0], 156, w, 22, "Air +\npsychrometrics", "rho, mu, k, c_p, Pr\nw, dew point\np_sat, rho_v,sat, h_fg", PALETTE["shared_a"], body_fs=7.4)
    box(ax, x_cols[1], 156, w, 22, "Bin loss +\nthermal network", "h_ext, U_wall\nUA_tot, UA_binf, UA_tinf, UA_bt\nC_total, Bi, Q_req", PALETTE["shared_b"], body_fs=7.0)
    box(ax, x_cols[2], 156, w, 22, "Branch partition\nvariables", "Vdot_tube, mdot_tube\nA_hole, u_hole\nreleased mdot\nsegment dx", PALETTE["shared_c"], body_fs=7.2)
    box(ax, x_cols[3], 156, w, 22, "Transport + loss\nclosures", "Re, Nu_wire, Nu_tube, f\nU_heater, U_aer\nK_contr, K_exp, dP_dist\nh_mass from Lewis", PALETTE["shared_d"], body_fs=6.8)
    box(ax, x_cols[4], 156, w, 22, "Electrical derived\nvariables", "L_wire, R_tube\nI_tube, P_tube\nq_gen,i", PALETTE["shared_e"], body_fs=7.4)

    # Heating branch
    box(ax, x_cols[0], 108, w, 24, "Heater tube\nmarch", "q_loss,i = U_heater P dx (T_air - T_inf)\nT_air,heater(x), T_wire(x)\ndP_heater", PALETTE["heat"], body_fs=6.9)
    box(ax, x_cols[1], 108, w, 24, "Heating aeration\nmarch", "q_wall,i + q_jet,i - q_evap,i\nT_air,aer(x), mdot_rem(x)\ndP_aer + dP_exp", PALETTE["heat"], body_fs=6.9)
    box(ax, x_cols[2], 108, w, 24, "Heating bed\nstate", "Q_bed,total -> T_b, T_t\nQ_b = f_b Q_bed\nQ_t = (1-f_b) Q_bed", PALETTE["heat"], body_fs=7.1)
    box(ax, x_cols[3], 108, w, 24, "Heating residuals\n+ duty", "I, T_wire,max, T_air,out\ndP, water, u_h\ndutyNeeded = Q_req / Q_bed,total", PALETTE["heat"], body_fs=6.9)
    box(ax, x_cols[4], 108, w, 24, "Reported heating\nreference", "Pareto / fallback point\nV*, flow*, Q_bed*\nP*, T_wire*, T_air,out*", PALETTE["heat"], body_fs=7.0)

    # Cooling branch
    box(ax, x_cols[0], 62, w, 24, "Spot cooler +\npsychrometrics", "Q_req,cool -> Q_actual,cool\nT_supply, m_cond, w_out\nP_spot", PALETTE["cool"], body_fs=7.1)
    box(ax, x_cols[1], 62, w, 24, "Assist blower +\nstatic limit", "V_free\ndP_avail(V)\nP_blower", PALETTE["cool"], body_fs=7.3)
    box(ax, x_cols[2], 62, w, 24, "Summer branch\nmarch", "q_wall,i + q_jet,i - q_evap,i\nQ_bed,total(<0), water, dP\nT_air,out", PALETTE["cool"], body_fs=6.9)
    box(ax, x_cols[3], 62, w, 24, "Cooling bed\nclosure", "R(T_bed,ref) = 0\nT_b, T_t, Q_cool\nshortfall + hard residuals", PALETTE["cool"], body_fs=7.0)
    box(ax, x_cols[4], 62, w, 24, "Reported cooling\nreference", "Pareto / fallback point\nflow*, Q_cool*\nP_spot* + P_blower*", PALETTE["cool"], body_fs=7.0)

    # Annual + post
    box(ax, x_cols[0], 18, w, 24, "Representative\nclimate year", "kernel baseline\nfreeze override\nheat-wave override\nT_inf(d)", PALETTE["annual"], body_fs=7.0)
    box(ax, x_cols[1], 18, w, 24, "Daily required\nloads", "Q_req,h(d) = UA bin max(T_h,set - T_inf, 0)\nQ_req,c(d) = UA bin max(T_inf - T_c,set, 0)", PALETTE["annual"], body_fs=6.5)
    box(ax, x_cols[2], 18, w, 24, "Daily heating\nre-solve", "Q_cap,h(d), P_h(d)\nwater_h(d)\nfrom heating ref + T_inf(d)", PALETTE["annual"], body_fs=7.0)
    box(ax, x_cols[3], 18, w, 24, "Daily cooling\nre-solve", "Q_cap,c(d), P_c(d)\nwater_c(d)\nfrom cooling ref + T_inf(d)", PALETTE["annual"], body_fs=7.0)
    box(ax, x_cols[4], 18, w, 24, "Duty + daily\ntotals", "D_h(d), D_c(d), Q_bed(d)\nT_b(d), T_t(d)\nkWh(d), cost(d), water(d), unmet(d)", PALETTE["annual"], body_fs=6.8)
    box(ax, 154, 3, 112, 8.0, "Post-processing appendix", "selected branch states -> plume radius/contact area -> habitat exclusion + cross section", PALETTE["post"], title_fs=9.2, body_fs=6.9)

    # Input to shared arrows
    arrow(ax, center_bottom(x_cols[0], 200, w, 24), center_top(x_cols[0], 156, w, 22), color=PALETTE["blue"])
    arrow(ax, center_bottom(x_cols[1], 200, w, 24), center_top(x_cols[1], 156, w, 22), color=PALETTE["green"])
    arrow(ax, center_bottom(x_cols[2], 200, w, 24), center_top(x_cols[2], 156, w, 22), color=PALETTE["orange"])
    arrow(ax, right_mid(x_cols[2], 200, w, 24), left_mid(x_cols[3], 156, w, 22), color=PALETTE["orange"], rad=-0.08)
    arrow(ax, center_bottom(x_cols[3], 200, w, 24), center_top(x_cols[4], 156, w, 22), color=PALETTE["red"], rad=-0.04)
    arrow(ax, center_bottom(x_cols[3], 200, w, 24), center_top(x_cols[3], 156, w, 22), color=PALETTE["red"])
    arrow(ax, center_bottom(x_cols[4], 200, w, 24), center_top(x_cols[0], 62, w, 24), color=PALETTE["teal"], rad=0.08)
    arrow(ax, center_bottom(x_cols[4], 200, w, 24), center_top(x_cols[1], 62, w, 24), color=PALETTE["teal"], rad=0.03)
    arrow(ax, center_bottom(x_cols[0], 200, w, 24), center_top(x_cols[0], 18, w, 24), color=PALETTE["gray"], rad=0.05)

    # Shared to heating
    arrow(ax, center_bottom(x_cols[0], 156, w, 22), center_top(x_cols[0], 108, w, 24), color=PALETTE["blue"])
    arrow(ax, right_mid(x_cols[0], 156, w, 22), left_mid(x_cols[1], 108, w, 24), color=PALETTE["blue"], rad=-0.06)
    arrow(ax, center_bottom(x_cols[2], 156, w, 22), center_top(x_cols[0], 108, w, 24), color=PALETTE["orange"], rad=0.06)
    arrow(ax, center_bottom(x_cols[2], 156, w, 22), center_top(x_cols[1], 108, w, 24), color=PALETTE["orange"])
    arrow(ax, center_bottom(x_cols[3], 156, w, 22), center_top(x_cols[0], 108, w, 24), color=PALETTE["teal"], rad=0.08)
    arrow(ax, center_bottom(x_cols[3], 156, w, 22), center_top(x_cols[1], 108, w, 24), color=PALETTE["teal"], rad=0.02)
    arrow(ax, center_bottom(x_cols[4], 156, w, 22), center_top(x_cols[0], 108, w, 24), color=PALETTE["red"])
    arrow(ax, center_bottom(x_cols[1], 156, w, 22), center_top(x_cols[2], 108, w, 24), color=PALETTE["green"])
    arrow(ax, center_bottom(x_cols[1], 156, w, 22), center_top(x_cols[3], 108, w, 24), color=PALETTE["green"], rad=-0.05)

    # Heating internal
    arrow(ax, right_mid(x_cols[0], 108, w, 24), left_mid(x_cols[1], 108, w, 24), color=PALETTE["red"])
    arrow(ax, right_mid(x_cols[1], 108, w, 24), left_mid(x_cols[2], 108, w, 24), color=PALETTE["red"])
    arrow(ax, right_mid(x_cols[2], 108, w, 24), left_mid(x_cols[3], 108, w, 24), color=PALETTE["red"])
    arrow(ax, right_mid(x_cols[3], 108, w, 24), left_mid(x_cols[4], 108, w, 24), color=PALETTE["red"])

    # Shared to cooling
    arrow(ax, center_bottom(x_cols[0], 156, w, 22), center_top(x_cols[0], 62, w, 24), color=PALETTE["blue"], rad=0.03)
    arrow(ax, center_bottom(x_cols[0], 156, w, 22), center_top(x_cols[2], 62, w, 24), color=PALETTE["blue"], rad=-0.10)
    arrow(ax, center_bottom(x_cols[2], 156, w, 22), center_top(x_cols[2], 62, w, 24), color=PALETTE["orange"])
    arrow(ax, center_bottom(x_cols[3], 156, w, 22), center_top(x_cols[2], 62, w, 24), color=PALETTE["teal"], rad=0.06)
    arrow(ax, center_bottom(x_cols[1], 156, w, 22), center_top(x_cols[3], 62, w, 24), color=PALETTE["green"], rad=-0.03)

    # Cooling internal
    arrow(ax, right_mid(x_cols[0], 62, w, 24), left_mid(x_cols[2], 62, w, 24), color=PALETTE["teal"], rad=-0.10)
    arrow(ax, right_mid(x_cols[1], 62, w, 24), left_mid(x_cols[3], 62, w, 24), color=PALETTE["teal"], rad=-0.06)
    arrow(ax, right_mid(x_cols[2], 62, w, 24), left_mid(x_cols[3], 62, w, 24), color=PALETTE["teal"])
    arrow(ax, right_mid(x_cols[3], 62, w, 24), left_mid(x_cols[4], 62, w, 24), color=PALETTE["teal"])

    # Shared/input to annual
    arrow(ax, center_bottom(x_cols[1], 156, w, 22), center_top(x_cols[1], 18, w, 24), color=PALETTE["green"], rad=0.08)
    arrow(ax, center_bottom(x_cols[0], 156, w, 22), center_top(x_cols[3], 18, w, 24), color=PALETTE["blue"], rad=-0.16)
    arrow(ax, center_bottom(x_cols[1], 156, w, 22), center_top(x_cols[4], 18, w, 24), color=PALETTE["green"], rad=-0.10)

    # References into annual and post-processing
    arrow(ax, center_bottom(x_cols[4], 108, w, 24), center_top(x_cols[2], 18, w, 24), color=PALETTE["red"], rad=0.12)
    arrow(ax, center_bottom(x_cols[4], 62, w, 24), center_top(x_cols[3], 18, w, 24), color=PALETTE["teal"], rad=-0.04)
    arrow(ax, center_bottom(x_cols[4], 108, w, 24), center_top(154, 3, 112, 8.0), color=PALETTE["gray"], rad=-0.02)
    arrow(ax, center_bottom(x_cols[4], 62, w, 24), center_top(154, 3, 112, 8.0), color=PALETTE["gray"], rad=0.02)

    # Annual internal
    arrow(ax, right_mid(x_cols[0], 18, w, 24), left_mid(x_cols[1], 18, w, 24), color=PALETTE["gold"])
    arrow(ax, right_mid(x_cols[1], 18, w, 24), left_mid(x_cols[4], 18, w, 24), color=PALETTE["gold"], rad=-0.18)
    arrow(ax, right_mid(x_cols[2], 18, w, 24), left_mid(x_cols[4], 18, w, 24), color=PALETTE["gold"], rad=-0.08)
    arrow(ax, right_mid(x_cols[3], 18, w, 24), left_mid(x_cols[4], 18, w, 24), color=PALETTE["gold"], rad=0.02)

    ax.text(
        135,
        0.9,
        "Arrows denote explicit variable dependence in the live model; grouped boxes collect tightly coupled equation sets rather than function-call boundaries.",
        ha="center",
        va="bottom",
        fontsize=8.2,
        color=PALETTE["muted"],
    )

    svg_path = OUT_DIR / "strict_physics_dependency_map_current.svg"
    png_path = OUT_DIR / "strict_physics_dependency_map_current.png"
    fig.savefig(svg_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=220, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved {svg_path}")
    print(f"Saved {png_path}")


if __name__ == "__main__":
    build_map()

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle


OUT_DIR = Path(r"C:\Users\cezar\CHEG4143W_HeatingModel\heating_diagrams")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def setup_fig(title, size=(14, 8)):
    fig, ax = plt.subplots(figsize=size)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 70)
    ax.axis("off")
    fig.suptitle(title, fontsize=18, fontweight="bold", y=0.98)
    return fig, ax


def label_box(ax, xy, w, h, text, fc, ec="#1f2937", fontsize=11, lw=2, rounded=True):
    if rounded:
        patch = FancyBboxPatch(
            xy,
            w,
            h,
            boxstyle="round,pad=0.6,rounding_size=1.4",
            linewidth=lw,
            edgecolor=ec,
            facecolor=fc,
        )
    else:
        patch = Rectangle(xy, w, h, linewidth=lw, edgecolor=ec, facecolor=fc)
    ax.add_patch(patch)
    ax.text(xy[0] + w / 2, xy[1] + h / 2, text, ha="center", va="center", fontsize=fontsize)
    return patch


def arrow(ax, p0, p1, text=None, color="#111827", lw=2.2, ms=16, text_offset=(0, 0)):
    arr = FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=ms, linewidth=lw, color=color)
    ax.add_patch(arr)
    if text:
        xm = 0.5 * (p0[0] + p1[0]) + text_offset[0]
        ym = 0.5 * (p0[1] + p1[1]) + text_offset[1]
        ax.text(xm, ym, text, fontsize=10, ha="center", va="center", color=color)
    return arr


def save(fig, name):
    png = OUT_DIR / f"{name}.png"
    svg = OUT_DIR / f"{name}.svg"
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)


def diagram_01_system_overview():
    fig, ax = setup_fig("01. Modeled System Overview")

    label_box(ax, (3, 6), 94, 58, "Greenhouse Winter Air Domain\nStagnant external air used for bin heat-loss model", "#eef7ff")
    ax.text(8, 59, "Outside design air used in model: greenhouse air at -15 C", fontsize=11, color="#1d4ed8")

    label_box(ax, (12, 15), 64, 36, "Vermicomposter Bin\n0.61 m x 0.61 m x 1.22 m", "#fef3c7")
    ax.text(16, 45, "Two-node bed model", fontsize=11, fontweight="bold")
    ax.text(16, 41, "Bottom node receives 60% of heater input", fontsize=10)
    ax.text(16, 38, "Top node receives 40% of heater input", fontsize=10)
    ax.text(16, 35, "Heat loss to greenhouse through insulated walls", fontsize=10)

    # Four parallel tubes
    yvals = [22, 28, 34, 40]
    for i, y in enumerate(yvals, start=1):
        label_box(ax, (22, y - 1.4), 46, 2.8, "", "#d1fae5", rounded=False, lw=1.5)
        for x in [30, 38, 46, 54, 62]:
            ax.add_patch(Circle((x, y), 0.35, color="#059669"))
        ax.text(70, y, f"Aeration tube {i}", va="center", fontsize=9)

    # Heater branches
    for i, y in enumerate(yvals, start=1):
        label_box(ax, (4, y - 2.0), 12, 4.0, f"Heater branch {i}\nNiCr coil", "#fecaca", fontsize=9)
        arrow(ax, (16, y), (22, y), "branch flow", color="#b91c1c", text_offset=(0, 1.5))

    label_box(ax, (3, 54), 18, 7, "Shared DC supply\nVoltage sweep", "#e9d5ff")
    for y in yvals:
        arrow(ax, (21, 57.5), (10, y + 1.5), color="#7c3aed", lw=1.6, ms=12)

    label_box(ax, (79, 18), 14, 24, "PID / Sensors\nBottom thermistor\nTop thermistor\nVoltage command", "#ddd6fe")
    arrow(ax, (76, 33), (79, 33), "T_bottom", color="#4f46e5", text_offset=(0, 2))
    arrow(ax, (76, 25), (79, 25), "T_top", color="#4f46e5", text_offset=(0, -2))
    arrow(ax, (86, 42), (21, 57.5), "V command", color="#4f46e5", text_offset=(0, 2))

    ax.text(
        50,
        10,
        "Total flow is split across parallel branches. Adding tubes improves distribution and reduces local branch outlet temperature.",
        ha="center",
        fontsize=10,
        color="#374151",
    )

    save(fig, "01_system_overview")


def diagram_02_heater_branch():
    fig, ax = setup_fig("02. Heater Branch Model: Joule Heating to Moving Air")

    label_box(ax, (4, 28), 12, 12, "Cold branch air\nfrom greenhouse\nT_in", "#dbeafe")
    label_box(ax, (20, 24), 56, 20, "Heater tube branch", "#fee2e2")
    label_box(ax, (80, 28), 14, 12, "Warm branch air\nto aeration tube\nT_out <= 30 C", "#dcfce7")

    arrow(ax, (16, 34), (20, 34), "m_dot_branch", color="#2563eb", text_offset=(0, 2))
    arrow(ax, (76, 34), (80, 34), "same branch flow", color="#2563eb", text_offset=(0, 2))

    # Wire helix impression
    xs = [24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]
    ys = [34, 38, 30, 38, 30, 38, 30, 38, 30, 38, 30, 38, 34]
    ax.plot(xs, ys, color="#b91c1c", linewidth=3)
    ax.text(49, 42, "16 AWG NiCr 80/20 helical wire", ha="center", fontsize=11, fontweight="bold")

    ax.text(49, 18, "Per-branch electrical model", ha="center", fontsize=11, fontweight="bold")
    ax.text(49, 14.5, "R_branch = rho(T) * L_wire / A_wire", ha="center", fontsize=11)
    ax.text(49, 11.5, "I_branch = V / R_branch    and    P_branch = V * I_branch", ha="center", fontsize=11)

    ax.text(49, 51, "Convective wire model used in code", ha="center", fontsize=11, fontweight="bold")
    ax.text(49, 47.5, "Re_wire from cylinder crossflow approximation on the helix", ha="center", fontsize=10)
    ax.text(49, 44.5, "Nu_wire = Churchill-Bernstein", ha="center", fontsize=10)
    ax.text(49, 41.5, "h_wire = Nu_wire * k_air / D_wire", ha="center", fontsize=10)

    ax.text(49, 58, "Axial air energy balance inside heater tube", ha="center", fontsize=11, fontweight="bold")
    ax.text(49, 54.5, "m_dot * c_p * dT = q_gen - q_loss", ha="center", fontsize=11)
    ax.text(49, 4.5, "The model also subtracts shell heat loss from the heater tube to the greenhouse air.", ha="center", fontsize=10)

    save(fig, "02_heater_branch")


def diagram_03_aeration_tube():
    fig, ax = setup_fig("03. Aeration Tube Model: Wall Heat + Released Jet Enthalpy")

    label_box(ax, (6, 28), 14, 12, "From heater branch\nT_out, m_dot_branch", "#dcfce7")
    label_box(ax, (24, 26), 58, 16, "Perforated aeration tube inside bed", "#d1fae5")
    arrow(ax, (20, 34), (24, 34), "abrupt expansion", color="#059669", text_offset=(0, 2))

    # Holes and release arrows
    for x in [34, 42, 50, 58, 66, 74]:
        ax.add_patch(Circle((x, 34), 0.35, color="#047857"))
        arrow(ax, (x, 34), (x, 24), color="#10b981", lw=1.5, ms=10)

    ax.text(53, 43, "Distributed release along tube", ha="center", fontsize=11, fontweight="bold")
    ax.text(53, 20, "Bed / substrate domain at T_bed", ha="center", fontsize=11, color="#14532d")

    ax.text(53, 56, "Per-segment thermal model", ha="center", fontsize=11, fontweight="bold")
    ax.text(53, 52.5, "q_wall = U * P * dx * (T_air - T_bed)", ha="center", fontsize=11)
    ax.text(53, 49.0, "q_jet = m_dot_release * c_p * (T_air - T_bed)", ha="center", fontsize=11)
    ax.text(53, 45.5, "Q_to_bed = sum(q_wall + q_jet) along the tube", ha="center", fontsize=11)

    ax.text(53, 10, "Each branch solves one aeration tube; total bed heating is multiplied by the number of parallel branches.", ha="center", fontsize=10)
    ax.text(53, 6.5, "Pressure drop includes heater friction, sudden expansion, aeration-tube friction, and a header minor loss.", ha="center", fontsize=10)

    save(fig, "03_aeration_tube_heat_release")


def diagram_04_bin_thermal_model():
    fig, ax = setup_fig("04. Bin Thermal Model: Lumped Heat Loss + Two-Node Bed")

    label_box(ax, (5, 14), 22, 42, "Greenhouse air\nT_inf\nNatural convection\noutside bin wall", "#eef2ff")
    label_box(ax, (35, 14), 22, 18, "Bottom bed node\nC_bottom\nT_bottom", "#fde68a")
    label_box(ax, (35, 38), 22, 18, "Top bed node\nC_top\nT_top", "#fde68a")
    label_box(ax, (67, 14), 22, 42, "Insulated wall path\nU*A to greenhouse\nfor each node", "#e0f2fe")

    arrow(ax, (27, 23), (35, 23), "Q_loss,bottom", color="#2563eb", text_offset=(0, 2))
    arrow(ax, (27, 47), (35, 47), "Q_loss,top", color="#2563eb", text_offset=(0, 2))
    arrow(ax, (57, 23), (67, 23), "same wall path", color="#0369a1", text_offset=(0, -2))
    arrow(ax, (57, 47), (67, 47), "same wall path", color="#0369a1", text_offset=(0, -2))
    arrow(ax, (46, 32), (46, 38), "UA_internal", color="#92400e", text_offset=(6, 0))
    arrow(ax, (46, 38), (46, 32), color="#92400e")

    arrow(ax, (89, 24), (57, 24), "heater input mostly to bottom", color="#b91c1c", text_offset=(0, 2))
    arrow(ax, (89, 46), (57, 46), "remainder to top", color="#b91c1c", text_offset=(0, 2))

    ax.text(50, 61, "Lumped external loss model used for sizing", ha="center", fontsize=11, fontweight="bold")
    ax.text(50, 57.5, "Q_loss = UA * (T_bin - T_inf)", ha="center", fontsize=11)
    ax.text(50, 54.0, "tau = C / UA", ha="center", fontsize=11)
    ax.text(50, 50.5, "T(t) = T_inf + (T0 - T_inf) exp(-t/tau)", ha="center", fontsize=11)

    ax.text(50, 7.5, "The code also evaluates Bi = h * Lc / k. In the verified run, Bi > 0.1, so lumped loss is acceptable for sizing but not for internal gradients.", ha="center", fontsize=10)

    save(fig, "04_bin_lumped_and_two_node_model")


def diagram_05_constraints():
    fig, ax = setup_fig("05. Main Design Constraint Identified by the Model")

    label_box(ax, (7, 18), 26, 34, "Winter requirement\nQ_loss ~ 61.6 W\nat T_bed = 20 C", "#fee2e2")
    label_box(ax, (37, 18), 26, 34, "Conservative worm-safety cap\nT_injected <= 30 C\nso DeltaT_use = 10 C", "#fef3c7")
    label_box(ax, (67, 18), 26, 34, "Thermodynamic consequence\nQ = m_dot c_p (T_injected - T_bed)", "#dcfce7")

    arrow(ax, (33, 35), (37, 35), color="#111827")
    arrow(ax, (63, 35), (67, 35), color="#111827")

    ax.text(50, 61, "Why four tubes helped but did not achieve 20% duty", ha="center", fontsize=12, fontweight="bold")
    ax.text(50, 56, "More tubes reduce local branch outlet temperature and improve spatial distribution.", ha="center", fontsize=11)
    ax.text(50, 52, "They do not reduce the minimum total airflow required at a fixed injected-air temperature cap.", ha="center", fontsize=11)

    ax.text(50, 11, "Verified model result", ha="center", fontsize=11, fontweight="bold")
    ax.text(50, 7.5, "Minimum total flow at 30 C cap: ~268 L/min at 100% duty, ~1342 L/min at 20% duty", ha="center", fontsize=11)

    save(fig, "05_core_constraint_and_duty_cycle")


def main():
    diagram_01_system_overview()
    diagram_02_heater_branch()
    diagram_03_aeration_tube()
    diagram_04_bin_thermal_model()
    diagram_05_constraints()
    print(f"Saved diagrams to {OUT_DIR}")


if __name__ == "__main__":
    main()

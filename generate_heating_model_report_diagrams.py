from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle


OUT_DIR = Path(r"C:\Users\cezar\CHEG4143W_HeatingModel\heating_report_diagrams_v2")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def fig_ax(title, size=(16, 10), xlim=(0, 100), ylim=(0, 100)):
    fig, ax = plt.subplots(figsize=size)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.axis("off")
    fig.suptitle(title, fontsize=22, fontweight="bold", y=0.985)
    return fig, ax


def box(ax, x, y, w, h, text, fc="#ffffff", ec="#1f2937", lw=2.2, fs=12, rounded=True):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.5,rounding_size=1.5" if rounded else "square,pad=0.0",
        linewidth=lw,
        edgecolor=ec,
        facecolor=fc,
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs)
    return patch


def line_arrow(ax, p0, p1, text=None, color="#111827", lw=2.4, ms=18, text_xy=None, fs=11):
    arr = FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=ms, linewidth=lw, color=color)
    ax.add_patch(arr)
    if text:
        tx, ty = text_xy if text_xy is not None else ((p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2)
        ax.text(tx, ty, text, fontsize=fs, ha="center", va="center", color=color)
    return arr


def note(ax, x, y, text, fs=11, color="#374151", ha="left"):
    ax.text(x, y, text, fontsize=fs, color=color, ha=ha, va="center")


def save(fig, stem):
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=240, bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def system_overview():
    fig, ax = fig_ax("01. Report Diagram: System Overview", size=(18, 10))
    box(ax, 2, 6, 96, 88, "", fc="#edf6ff", fs=16)
    note(ax, 8, 88, "External boundary used in the heat-loss model: stagnant greenhouse air, design temperature = -15 C", fs=13, color="#1d4ed8")

    box(ax, 20, 18, 56, 54, "", fc="#fbf3cf", fs=18)
    note(ax, 48, 64, "Vermicomposter bin", fs=18, ha="center")
    note(ax, 48, 59.5, "Modeled internally as a bottom node and a top node with conductive coupling.", fs=12.5, ha="center")
    note(ax, 48, 55.0, "Heater input enters through four parallel aeration branches.", fs=12.5, ha="center")
    note(ax, 48, 51.8, "Heat leaves through the insulated wall to greenhouse air.", fs=12.5, ha="center")

    yvals = [24, 32, 40, 48]
    for i, y in enumerate(yvals, start=1):
        box(ax, 21, y - 1.8, 46, 3.6, "", fc="#cdeee0", ec="#0f766e", lw=1.6, rounded=False)
        for x in [30, 38, 46, 54, 62]:
            ax.add_patch(Circle((x, y), 0.35, facecolor="#059669", edgecolor="#047857"))
        note(ax, 69, y, f"Aeration tube {i}", fs=11)

    box(ax, 4, 72, 14, 10, "Shared DC supply", fc="#e9ddff", fs=13)
    box(ax, 4, 20, 12, 40, "Four parallel\nheater branches\n(NiCr coil in each)", fc="#f8d1d1", fs=13)
    for y in yvals:
        line_arrow(ax, (16, y), (21, y), text="branch flow", color="#b91c1c", text_xy=(18.5, y + 2.4), fs=11)
        line_arrow(ax, (11, 72), (11, y + 2.0), color="#6d28d9", lw=1.8, ms=12)

    box(ax, 80, 28, 14, 24, "PID / sensors\nbottom thermistor\ntop thermistor\nvoltage command", fc="#dcd6f7", fs=14)
    line_arrow(ax, (76, 31), (80, 31), text="T_bottom", color="#4f46e5", text_xy=(77.5, 34.5), fs=11)
    line_arrow(ax, (76, 45), (80, 45), text="T_top", color="#4f46e5", text_xy=(77.5, 48.5), fs=11)
    line_arrow(ax, (87, 52), (87, 74), color="#4f46e5")
    line_arrow(ax, (87, 74), (18, 77), text="voltage command", color="#4f46e5", text_xy=(56, 79), fs=12)

    note(ax, 50, 12, "Total system flow is split across four equal branches in the current model.", fs=13, ha="center")
    save(fig, "01_report_system_overview")


def heater_branch():
    fig, ax = fig_ax("02. Report Diagram: One Heater Branch", size=(18, 10))
    box(ax, 4, 38, 14, 12, "Cold branch air\nfrom greenhouse", fc="#dbeafe", fs=14)
    box(ax, 24, 28, 52, 22, "Heater tube branch", fc="#f8d7d7", fs=18)
    box(ax, 82, 38, 14, 12, "Warm branch air\nto aeration tube", fc="#d9f2df", fs=14)
    line_arrow(ax, (18, 44), (24, 44), text="m_dot,branch", color="#2563eb", text_xy=(21, 47), fs=12)
    line_arrow(ax, (76, 44), (82, 44), text="T_out,branch", color="#2563eb", text_xy=(79, 47), fs=12)

    xs = [28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]
    ys = [44, 50, 38, 50, 38, 50, 38, 50, 38, 50, 38, 44]
    ax.plot(xs, ys, color="#c02626", linewidth=3.0)
    note(ax, 50, 53, "16 AWG NiCr 80/20 helix", fs=14, ha="center")

    box(ax, 6, 66, 28, 18, "Electrical submodel\n\nR_branch = rho(T) L / A\nI_branch = V / R_branch\nP_branch = V I_branch", fc="#f3e8ff", fs=14)
    box(ax, 36, 66, 28, 18, "Wire convection submodel\n\nRe_wire from crossflow on helix\nNu_wire = Churchill-Bernstein\nh_wire = Nu_wire k / D_wire", fc="#fee2e2", fs=13)
    box(ax, 66, 66, 28, 18, "Axial air energy balance\n\nm_dot c_p dT = q_gen - q_loss\nq_loss = U P dx (T_air - T_env)", fc="#e0f2fe", fs=13)

    line_arrow(ax, (20, 66), (33, 56), color="#7c3aed")
    line_arrow(ax, (50, 66), (50, 56), color="#b91c1c")
    line_arrow(ax, (80, 66), (67, 56), color="#0369a1")

    box(ax, 22, 8, 56, 12, "Outputs from one branch\n\nT_air,out,branch   T_wire,max   h_wire,mean   DeltaP_heater", fc="#ecfccb", fs=14)
    line_arrow(ax, (50, 28), (50, 20), color="#166534")
    save(fig, "02_report_heater_branch")


def aeration_branch():
    fig, ax = fig_ax("03. Report Diagram: One Aeration Branch", size=(18, 10))
    box(ax, 4, 42, 16, 12, "Input from heater\nT_air,in,branch\nm_dot,branch", fc="#dcfce7", fs=14)
    box(ax, 24, 34, 56, 16, "Perforated aeration tube in the bed", fc="#d1fae5", fs=18)
    line_arrow(ax, (20, 48), (24, 48), text="abrupt expansion", color="#059669", text_xy=(22, 52), fs=11)

    for x in [34, 42, 50, 58, 66, 74]:
        ax.add_patch(Circle((x, 42), 0.38, facecolor="#047857", edgecolor="#065f46"))
        line_arrow(ax, (x, 42), (x, 30), color="#10b981", lw=1.6, ms=10)

    box(ax, 10, 66, 24, 18, "Internal tube convection\n\nRe_tube\nNu_tube\nh_inside = Nu k / D", fc="#e0f2fe", fs=13)
    box(ax, 38, 66, 24, 18, "Wall transfer to bed\n\nU from inside h, wall k,\noutside h_bed\nq_wall = U P dx (T_air - T_bed)", fc="#fef3c7", fs=13)
    box(ax, 66, 66, 24, 18, "Released jet enthalpy\n\nq_jet = m_dot_release c_p\n(T_air - T_bed)", fc="#fee2e2", fs=13)

    line_arrow(ax, (22, 66), (34, 54), color="#0284c7")
    line_arrow(ax, (50, 66), (50, 54), color="#a16207")
    line_arrow(ax, (78, 66), (66, 54), color="#b91c1c")

    box(ax, 20, 10, 60, 14, "Branch outputs\n\nQ_to_bed,branch = sum(q_wall + q_jet)\nDeltaP_branch = DeltaP_heater + DeltaP_expansion + DeltaP_aeration", fc="#ecfccb", fs=14)
    line_arrow(ax, (52, 34), (50, 24), color="#166534")
    note(ax, 50, 4, "Total bed heating in the system is the branch value multiplied by the number of parallel branches.", fs=12, ha="center")
    save(fig, "03_report_aeration_branch")


def bin_model():
    fig, ax = fig_ax("04. Report Diagram: Bin Thermal Model", size=(18, 10))
    box(ax, 4, 22, 20, 46, "Greenhouse air\n\nT_inf\nnatural convection\noutside the insulated wall", fc="#e8eef8", fs=15)
    box(ax, 38, 48, 22, 18, "Top bed node\n\nC_top\nT_top", fc="#f8e28b", fs=15)
    box(ax, 38, 22, 22, 18, "Bottom bed node\n\nC_bottom\nT_bottom", fc="#f8e28b", fs=15)
    box(ax, 74, 22, 20, 46, "Heater input from\naeration network\n\n60% to bottom\n40% to top", fc="#fde2e2", fs=15)

    line_arrow(ax, (24, 57), (38, 57), text="Q_loss,top", color="#2563eb", text_xy=(31, 61), fs=12)
    line_arrow(ax, (24, 31), (38, 31), text="Q_loss,bottom", color="#2563eb", text_xy=(31, 35), fs=12)
    line_arrow(ax, (60, 57), (74, 57), text="heater split", color="#b91c1c", text_xy=(67, 61), fs=12)
    line_arrow(ax, (60, 31), (74, 31), text="heater split", color="#b91c1c", text_xy=(67, 35), fs=12)
    line_arrow(ax, (49, 40), (49, 48), text="UA_internal", color="#92400e", text_xy=(58, 44), fs=12)
    line_arrow(ax, (49, 48), (49, 40), color="#92400e")

    box(ax, 10, 76, 80, 14, "External loss sizing equations\n\nQ_loss = U A (T_bin - T_inf)      tau = C / (U A)      T(t) = T_inf + (T0 - T_inf) exp(-t/tau)", fc="#f0fdf4", fs=15)
    note(ax, 50, 10, "The code also checks Bi = h Lc / k. In the verified run, Bi > 0.1, so the lumped model is used only for loss sizing, not internal gradients.", fs=12, ha="center")
    save(fig, "04_report_bin_thermal_model")


def constraints():
    fig, ax = fig_ax("05. Report Diagram: Main Design Constraint", size=(18, 10))
    box(ax, 6, 28, 24, 30, "Winter loss requirement\n\nQ_loss ~ 61.6 W\nfor T_bed = 20 C", fc="#fee2e2", fs=16)
    box(ax, 38, 28, 24, 30, "Conservative safety cap\n\nT_injected <= 30 C\nso usable DeltaT = 10 C", fc="#fef3c7", fs=16)
    box(ax, 70, 28, 24, 30, "Thermodynamic limit\n\nQ = m_dot c_p\n(T_injected - T_bed)", fc="#dcfce7", fs=16)
    line_arrow(ax, (30, 43), (38, 43), color="#111827")
    line_arrow(ax, (62, 43), (70, 43), color="#111827")

    box(ax, 18, 66, 64, 18, "Verified model consequence\n\nMinimum total flow for 30 C injection cap:\n~268 L/min at 100% duty, ~1342 L/min at 20% duty", fc="#e0f2fe", fs=16)
    box(ax, 18, 8, 64, 12, "Adding more parallel tubes reduces local branch outlet temperature and improves heat distribution, but it does not remove the m_dot c_p DeltaT requirement.", fc="#f3f4f6", fs=14)
    save(fig, "05_report_main_constraint")


def equation_flow():
    fig, ax = fig_ax("06. Equation Dependency Map for the MATLAB Model", size=(30, 18), xlim=(0, 165), ylim=(0, 100))

    # Column 1: inputs
    box(ax, 3, 78, 20, 16, "Inputs\n\ngeometry\nmaterial properties\nvoltage sweep\nflow sweep\ntemperature limits", fc="#eef2ff", fs=15)
    box(ax, 3, 56, 20, 16, "Environment\n\nT_greenhouse\np_atm\nreference length", fc="#eef2ff", fs=15)
    box(ax, 3, 34, 20, 16, "Control settings\n\nKp Ki Kd\ntime step\ninitial T_bottom, T_top", fc="#eef2ff", fs=15)

    # Column 2: properties and base models
    box(ax, 28, 78, 22, 16, "Air properties\n\nmu(T): Sutherland\nk = mu c_p / Pr\nrho = p / (R T)", fc="#e0f2fe", fs=15)
    box(ax, 28, 56, 22, 16, "Natural convection outside bin\n\nRa = g beta DeltaT L^3 / (nu alpha)\nNu = Churchill-Chu\nh_ext = Nu k / L", fc="#e0f2fe", fs=14)
    box(ax, 28, 34, 22, 16, "Wire geometry and resistivity\n\nL_turn = sqrt((pi D)^2 + pitch^2)\nL_wire = N_turns L_turn\nR' = rho(T) / A", fc="#f3e8ff", fs=14)

    # Column 3: bin losses and heater branch
    box(ax, 55, 78, 22, 16, "Bin loss model\n\nU_wall = 1 / sum(R_i)\nQ_loss = U A (T_bin - T_inf)\ntau = C / (U A)\nBi = h Lc / k", fc="#fef3c7", fs=14)
    box(ax, 55, 56, 22, 16, "Branch electrical model\n\nR_branch = R' L_wire\nI_branch = V / R_branch\nP_branch = V I_branch", fc="#fee2e2", fs=15)
    box(ax, 55, 34, 22, 16, "Branch flow split\n\nm_dot_branch = m_dot_total / N_tubes\nu = m_dot / (rho A)", fc="#dcfce7", fs=15)

    # Column 4: heater tube physics
    box(ax, 82, 78, 22, 16, "Wire convection in heater tube\n\nRe_wire\nNu_wire = Churchill-Bernstein\nh_wire = Nu k / D_wire", fc="#fee2e2", fs=14)
    box(ax, 82, 56, 22, 16, "Internal tube convection\n\nRe_tube\nNu_tube = Hausen or Gnielinski\nh_inside = Nu k / D", fc="#fee2e2", fs=14)
    box(ax, 82, 34, 22, 16, "Heater energy balance\n\nq_loss = U P dx (T_air - T_env)\nq_gen from branch power\nm_dot c_p dT = q_gen - q_loss\nDeltaT_wire = q_gen / (h A)", fc="#fee2e2", fs=13)

    # Column 5: aeration + dp + aggregation
    box(ax, 109, 78, 22, 16, "Aeration tube heat transfer\n\nU from inside h, wall k, outside h_bed\nq_wall = U P dx (T_air - T_bed)\nq_jet = m_dot_release c_p (T_air - T_bed)", fc="#dcfce7", fs=13)
    box(ax, 109, 56, 22, 16, "Pressure drop model\n\nf = 64/Re or Churchill\nDeltaP_heater\nDeltaP_expansion\nDeltaP_aeration\nDeltaP_header", fc="#dcfce7", fs=13)
    box(ax, 109, 34, 22, 16, "System totals\n\nQ_bed,total = N_tubes Q_bed,branch\nP_total = N_tubes P_branch\nI_total = N_tubes I_branch\nR_eq = R_branch / N_tubes", fc="#dcfce7", fs=13)

    # Column 6: design selection + dynamic simulation
    box(ax, 136, 78, 22, 16, "Steady bed temperatures\n\nTwo-node linear balance\nbottom node + top node\nwith UA_internal coupling", fc="#fde68a", fs=14)
    box(ax, 136, 56, 22, 16, "Constraint check\n\nduty = Q_loss / Q_bed,total\nT_air,out <= limit\nT_wire,max <= limit\nDeltaP <= limit", fc="#fde68a", fs=14)
    box(ax, 136, 34, 22, 16, "PID transient\n\nerror = T_set - min(T_bottom, T_top)\nV_cmd from PID\nre-evaluate operating point\nupdate T_bottom and T_top", fc="#fde68a", fs=13)

    # Arrows between columns
    arrows = [
        ((23, 86), (28, 86)), ((23, 64), (28, 64)), ((23, 42), (28, 42)),
        ((50, 86), (55, 86)), ((50, 64), (55, 64)), ((50, 42), (55, 42)),
        ((77, 64), (82, 64)), ((77, 42), (82, 42)), ((77, 86), (82, 86)),
        ((104, 86), (109, 86)), ((104, 64), (109, 64)), ((104, 42), (109, 42)),
        ((131, 86), (136, 86)), ((131, 64), (136, 64)), ((131, 42), (136, 42)),
    ]
    for p0, p1 in arrows:
        line_arrow(ax, p0, p1, color="#111827", lw=2.0, ms=14)

    # Cross-links
    line_arrow(ax, (39, 78), (66, 72), text="h_ext feeds U_wall", color="#2563eb", text_xy=(52, 76), fs=11)
    line_arrow(ax, (39, 42), (66, 56), text="R' and L_wire feed branch R", color="#7c3aed", text_xy=(52, 50), fs=11)
    line_arrow(ax, (66, 34), (93, 78), text="branch velocity feeds Re_wire", color="#b91c1c", text_xy=(79, 58), fs=11)
    line_arrow(ax, (66, 34), (93, 56), text="branch velocity feeds Re_tube", color="#b91c1c", text_xy=(79, 46), fs=11)
    line_arrow(ax, (93, 56), (120, 78), text="inside h feeds aeration U", color="#0369a1", text_xy=(106, 70), fs=11)
    line_arrow(ax, (93, 34), (120, 34), text="T_air,out and branch heat", color="#166534", text_xy=(106, 30), fs=11)
    line_arrow(ax, (66, 86), (147, 56), text="Q_loss feeds duty check", color="#a16207", text_xy=(107, 78), fs=11)
    line_arrow(ax, (120, 34), (147, 86), text="Q_bed,total feeds bed temperatures", color="#059669", text_xy=(134, 60), fs=11)
    line_arrow(ax, (147, 86), (147, 34), text="temperatures feed PID", color="#92400e", text_xy=(154, 60), fs=11)

    note(ax, 82.5, 18, "Read left to right: inputs -> properties -> losses and branch models -> heater correlations -> aeration and totals -> constraints and PID.", fs=14, ha="center")
    note(ax, 82.5, 12, "This map matches the current MATLAB implementation and includes both physical equations and the computational sequence used in the design sweep.", fs=13, ha="center")

    save(fig, "06_equation_dependency_map")


def main():
    system_overview()
    heater_branch()
    aeration_branch()
    bin_model()
    constraints()
    equation_flow()
    print(f"Saved improved report diagrams to {OUT_DIR}")


if __name__ == "__main__":
    main()

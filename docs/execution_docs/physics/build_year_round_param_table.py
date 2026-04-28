from __future__ import annotations

import re
from collections import OrderedDict
from pathlib import Path


ROOT = Path(__file__).resolve().parent
WALKTHROUGH = ROOT / "section18_report_walkthrough.tex"
CONNECTION_GUIDE = ROOT / "section18_variable_connection_guide.tex"
OUT_MD = ROOT / "year_round_operation_parameter_table.md"


def clean_symbol(sym: str) -> str:
    s = sym.strip()
    s = s.replace(r"\), \(", ", ")
    s = s.replace(r"\),\(", ", ")
    s = s.replace(r"$", "")
    s = re.sub(r"\s*,\s*", ",", s)
    s = re.sub(r"\s+", " ", s)
    return s.strip(" ,")


def clean_definition(defn: str) -> str:
    d = defn.strip()
    d = re.sub(r"\\texttt\{([^}]*)\}", r"`\1`", d)
    d = d.replace(r"\(", "$").replace(r"\)", "$")
    d = re.sub(r"\s+", " ", d)
    return d


def extract_walkthrough_items(text: str) -> list[tuple[str, str]]:
    # Captures lines like: \item \(symbol\): definition...
    pattern = re.compile(r"\\item\s+\\\((.+?)\\\):\s*([^\n]+)")
    rows: list[tuple[str, str]] = []
    for m in pattern.finditer(text):
        sym = clean_symbol(m.group(1))
        defn = clean_definition(m.group(2))
        if sym and defn:
            rows.append((sym, defn))
    return rows


def extract_connection_rows(text: str) -> list[tuple[str, str]]:
    # Captures table lines like: $symbol$ & role & ...
    rows: list[tuple[str, str]] = []
    for line in text.splitlines():
        m = re.match(r"^\$(.+?)\$\s*&\s*(.+?)\s*&", line.strip())
        if not m:
            continue
        sym = clean_symbol(m.group(1))
        role = clean_definition(m.group(2))
        if "&" in sym:
            continue
        if sym and role:
            rows.append((sym, role))
    return rows


def add_or_keep(d: OrderedDict[str, str], sym: str, defn: str) -> None:
    if sym not in d:
        d[sym] = defn


def main() -> None:
    walkthrough_text = WALKTHROUGH.read_text(encoding="utf-8", errors="ignore")
    connection_text = CONNECTION_GUIDE.read_text(encoding="utf-8", errors="ignore")

    entries: OrderedDict[str, str] = OrderedDict()

    for sym, defn in extract_walkthrough_items(walkthrough_text):
        add_or_keep(entries, sym, defn)
    for sym, defn in extract_connection_rows(connection_text):
        add_or_keep(entries, sym, defn)

    # Explicit aliases and symbols requested by user and commonly missed by pattern extraction.
    manual: list[tuple[str, str]] = [
        (r"s_{h,m,q}", "Feasible heat-wave spell start day for spell $q$ in month $m$ after non-overlap and month-bound clipping."),
        (r"\ell_{h,m,q}", "Heat-wave spell length (days) for spell $q$ in month $m$."),
        (r"\mathcal{B}_{h,m,q}", "Set of day-of-month indices occupied by heat-wave spell $q$ in month $m$."),
        (r"\mathcal{B}_{h,m}^{\mathrm{all}}", r"Union of all heat-wave day blocks in month $m$: $\bigcup_q \mathcal{B}_{h,m,q}$."),
        (r"d_{h,m,q}^{\mathrm{raw}}", "Raw ordered start day for heat-wave spell $q$ before feasibility projection."),
        (r"e_{h,m,q}", "Earliest feasible start day for heat-wave spell $q$."),
        (r"u_{h,m,q}", "Latest feasible start day for heat-wave spell $q$."),
        (r"\delta^h_{m,q}", "Alias for heat-wave spell anomaly magnitude; identical to $\Delta T_{h,m,q}$."),
        (r"\delta^f_{m,q}", "Alias for freeze-event anomaly magnitude; identical to $\Delta T_{f,m,q}$."),
        (r"\Delta T_{h,m,q}", "Allocated heat-wave anomaly magnitude for spell $q$ in month $m$ (same as $M_{h,m,q}$)."),
        (r"\Delta T_{f,m,q}", "Allocated freeze anomaly magnitude for event $q$ in month $m$."),
        (r"D(d)", "Active duty fraction for day $d$: equals $D_h(d)$ on heating days or $D_c(d)$ on cooling days."),
    ]
    for sym, defn in manual:
        add_or_keep(entries, clean_symbol(sym), clean_definition(defn))

    overrides: dict[str, str] = {
        clean_symbol(r"y_i"): "Reference monthly value used in the active kernel pass; either $\\bar T_i$ (baseline pass) or $\\sigma_i$ (spread pass).",
        clean_symbol(r"T_{h,\mathrm{set}}, T_{c,\mathrm{set}}"): "Resolved heating and cooling setpoints used by the year-round solver.",
        clean_symbol(r"\varnothing"): "Empty set. Also used with set operators: $A\\setminus B$ (difference), $A\\cup\\{x\\}$ (union/add), and $|A|$ (set cardinality).",
    }
    for sym, defn in overrides.items():
        entries[sym] = defn

    sorted_entries = sorted(entries.items(), key=lambda kv: kv[0].lower())

    lines: list[str] = []
    lines.append("# Year-round Operation Parameter Table")
    lines.append("")
    lines.append(
        "Definitions for symbols used across the full report Year-round Operation section "
        "(Sections 18.1 through 18.6), including explicit aliases such as "
        "$\\delta^h_{m,q}$ and $\\delta^f_{m,q}$."
    )
    lines.append("")
    lines.append("| Symbol | Definition |")
    lines.append("|---|---|")

    for sym, defn in sorted_entries:
        sym_cell = f"${sym}$" if "$" not in sym else sym
        defn_cell = defn.replace("|", r"\|")
        lines.append(f"| {sym_cell} | {defn_cell} |")

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {OUT_MD} with {len(sorted_entries)} rows.")


if __name__ == "__main__":
    main()

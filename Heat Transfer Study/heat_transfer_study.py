from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import multiprocessing as mp
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from scipy.interpolate import interp1d


ROOT_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = ROOT_DIR / "study_config.json"
HEATING_MODE = "heating"
COOLING_MODE = "summer_spot_cooler"
DEFAULT_PRIORITY_ORDER = (
    "current_excess_A",
    "wire_temp_excess_C",
    "air_outlet_excess_C",
    "pressure_drop_excess_Pa",
    "water_loss_excess_kg_day",
    "hole_velocity_excess_m_s",
    "bottom_temp_shortfall_C",
    "top_temp_shortfall_C",
    "temp_spread_excess_C",
    "minus_min_bed_temp_C",
    "power_W",
    "temp_spread_C",
    "minus_mean_bed_temp_C",
    "heat_shortfall",
)
HARD_CONSTRAINT_KEYS = (
    "current_excess_A",
    "wire_temp_excess_C",
    "air_outlet_excess_C",
    "pressure_drop_excess_Pa",
    "water_loss_excess_kg_day",
    "hole_velocity_excess_m_s",
    "bottom_temp_shortfall_C",
    "top_temp_shortfall_C",
    "temp_spread_excess_C",
)
DEFAULT_COOLING_PRIORITY_ORDER = (
    "assist_blower_pressure_excess_Pa",
    "pressure_drop_excess_Pa",
    "water_loss_excess_kg_day",
    "hole_velocity_excess_m_s",
    "bottom_temp_excess_C",
    "top_temp_excess_C",
    "temp_spread_excess_C",
    "max_bed_temp_C",
    "mean_bed_temp_C",
    "combined_cooling_power_W",
    "totalFlow_Lpm",
    "temp_spread_C",
    "minus_cooling_capacity_W",
)
PARTICLE_SURFACE_KEYS = ("gt40", "mm20to40", "mm10to20", "lt10")
PARTICLE_SURFACE_LABELS = {
    "gt40": "> 40 mm",
    "mm20to40": "20-40 mm",
    "mm10to20": "10-20 mm",
    "lt10": "< 10 mm",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Standalone heat-transfer study plotter. Reads exported model data "
            "and plotting/target instructions from JSON files and regenerates "
            "the study figures without calling MATLAB."
        )
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="Path to the study configuration JSON file.",
    )
    return parser.parse_args()


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8-sig") as fh:
        return json.load(fh)


def operating_mode(config: dict) -> str:
    mode = str(config.get("operating_mode", HEATING_MODE)).strip().lower()
    if mode == "summer_evaporative_cooling":
        return COOLING_MODE
    return mode or HEATING_MODE


def enabled_run_modes(config: dict) -> dict[str, bool]:
    run_modes = config.get("run_modes", {})
    if isinstance(run_modes, dict) and run_modes:
        return {
            "heating": bool(run_modes.get("heating", False)),
            "cooling": bool(run_modes.get("cooling", False)),
        }

    mode = operating_mode(config)
    return {
        "heating": mode == HEATING_MODE,
        "cooling": mode == COOLING_MODE,
    }


def python_parallel_config(config: dict) -> dict:
    parallel = dict(config.get("python_parallel", {}))
    parallel.setdefault("enabled", False)
    parallel.setdefault("workers", 0)
    parallel.setdefault("coolingSweep", True)
    parallel.setdefault("yearRound", True)
    return parallel


def model_overrides_signature(overrides: dict) -> str:
    encoded = json.dumps(overrides, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def resolve_path(root: Path, raw_path: str) -> Path:
    path = Path(raw_path)
    if path.is_absolute():
        return path
    return (root / path).resolve()


def as_float_array(values) -> np.ndarray:
    return np.array([np.nan if value is None else float(value) for value in values], dtype=float)


def sweep_array_from_config(sweep_cfg: dict, prefix: str, fallback_key: str) -> list[float] | None:
    direct = sweep_cfg.get(fallback_key)
    if isinstance(direct, list) and direct:
        return direct
    min_key = f"{prefix}_min"
    max_key = f"{prefix}_max"
    count_key = f"{prefix}_count"
    if min_key in sweep_cfg and max_key in sweep_cfg and count_key in sweep_cfg:
        count = max(2, int(sweep_cfg[count_key]))
        return np.linspace(float(sweep_cfg[min_key]), float(sweep_cfg[max_key]), count).tolist()
    return None


def strip_comment_keys(value):
    if isinstance(value, dict):
        cleaned = {}
        for key, item in value.items():
            if str(key).startswith("_comment"):
                continue
            cleaned[key] = strip_comment_keys(item)
        return cleaned
    if isinstance(value, list):
        return [strip_comment_keys(item) for item in value]
    return value


def apply_override_aliases(overrides: dict) -> dict:
    normalized = json.loads(json.dumps(overrides))

    bin_wall_alias = normalized.pop("bin_wall", None)
    if isinstance(bin_wall_alias, dict):
        bin_wall = dict(normalized.get("binWall", {}))
        if "sheetThickness_mm" in bin_wall_alias:
            bin_wall["sheetThickness_m"] = float(bin_wall_alias["sheetThickness_mm"]) / 1000.0
        if "insulationThickness_mm" in bin_wall_alias:
            bin_wall["insulationThickness_m"] = float(bin_wall_alias["insulationThickness_mm"]) / 1000.0
        for key in ("sheetK_W_mK", "insulationK_W_mK", "internalH_W_m2K"):
            if key in bin_wall_alias:
                bin_wall[key] = bin_wall_alias[key]
        normalized["binWall"] = bin_wall

    heater_alias = normalized.pop("heater_tube", None)
    if isinstance(heater_alias, dict):
        heater = dict(normalized.get("heaterTube", {}))
        if "length_m" in heater_alias:
            heater["length_m"] = float(heater_alias["length_m"])
        if "ID_mm" in heater_alias:
            heater["ID_m"] = float(heater_alias["ID_mm"]) / 1000.0
        if "OD_mm" in heater_alias:
            heater["OD_m"] = float(heater_alias["OD_mm"]) / 1000.0
        if "insulationThickness_mm" in heater_alias:
            heater["insulationThickness_m"] = float(heater_alias["insulationThickness_mm"]) / 1000.0
        for key in ("wallK_W_mK", "insulationK_W_mK", "externalH_W_m2K", "coilMeanD_m", "coilSpan_m"):
            if key in heater_alias:
                heater[key] = heater_alias[key]
        normalized["heaterTube"] = heater

    aeration = normalized.get("aeration")
    if isinstance(aeration, dict):
        aeration = dict(aeration)
        if "length_m" in aeration:
            aeration["length_m"] = float(aeration["length_m"])
        if "ID_mm" in aeration:
            aeration["ID_m"] = float(aeration["ID_mm"]) / 1000.0
        if "OD_mm" in aeration:
            aeration["OD_m"] = float(aeration["OD_mm"]) / 1000.0
        if "headerID_mm" in aeration:
            aeration["headerID_m"] = float(aeration["headerID_mm"]) / 1000.0
        if "splitterInletID_mm" in aeration:
            aeration["splitterInletID_m"] = float(aeration["splitterInletID_mm"]) / 1000.0
        if "branchConnectorID_mm" in aeration:
            aeration["branchConnectorID_m"] = float(aeration["branchConnectorID_mm"]) / 1000.0
        normalized["aeration"] = aeration

    return normalized


def normalized_model_overrides(config: dict) -> dict:
    overrides = apply_override_aliases(strip_comment_keys(json.loads(json.dumps(config.get("model_overrides", {})))))
    sweep_cfg = overrides.get("sweep")
    if isinstance(sweep_cfg, dict):
        voltage = sweep_array_from_config(sweep_cfg, "voltage_V", "voltage_V")
        flow = sweep_array_from_config(sweep_cfg, "totalFlow_Lpm", "totalFlow_Lpm")
        wire_d = sweep_array_from_config(sweep_cfg, "wireDiameter_mm", "wireDiameter_mm")
        coil_d = sweep_array_from_config(sweep_cfg, "coilMeanD_mm", "coilMeanD_mm")
        wire_len = sweep_array_from_config(sweep_cfg, "wireLength_m", "wireLength_m")
        if voltage is not None:
            sweep_cfg["voltage_V"] = voltage
        if flow is not None:
            sweep_cfg["totalFlow_Lpm"] = flow
        if wire_d is not None:
            sweep_cfg["wireDiameter_mm"] = wire_d
        if coil_d is not None:
            sweep_cfg["coilMeanD_mm"] = coil_d
        if wire_len is not None:
            sweep_cfg["wireLength_m"] = wire_len
    return overrides


def can_switch_active_bin_from_export(export_overrides: dict | None, current_overrides: dict) -> bool:
    if not isinstance(export_overrides, dict):
        return False

    def strip_bin_mode(overrides: dict) -> dict:
        cleaned = strip_comment_keys(json.loads(json.dumps(overrides)))
        bin_cfg = cleaned.get("bin")
        if isinstance(bin_cfg, dict):
            for key in ("openTop", "ventHoleCount", "ventHoleDiameter_m"):
                bin_cfg.pop(key, None)
            if not bin_cfg:
                cleaned.pop("bin", None)
        return cleaned

    return strip_bin_mode(export_overrides) == strip_bin_mode(current_overrides)


def auto_refresh_export_if_needed(config_path: Path, config: dict) -> None:
    auto_cfg = dict(config.get("auto_refresh_export", {}))
    if not auto_cfg.get("enabled", False):
        return

    data_json = resolve_path(config_path.parent, config["data_json"])
    current_overrides = normalized_model_overrides(config)
    current_signature = model_overrides_signature(current_overrides)

    export_signature = None
    export_snapshot = None
    payload = None
    if data_json.exists():
        try:
            payload = load_json(data_json)
            export_signature = payload.get("meta", {}).get("model_overrides_signature")
            export_snapshot = payload.get("meta", {}).get("model_overrides_snapshot")
        except Exception:
            export_signature = None

    if export_signature == current_signature:
        return

    if (
        auto_cfg.get("allow_bin_only_switch_from_export", True)
        and payload is not None
        and can_switch_active_bin_from_export(export_snapshot, current_overrides)
    ):
        desired_label = requested_active_configuration_label(config)
        if desired_label and any(
            same_configuration_label(desired_label, item.get("label", ""))
            for item in payload.get("configuration_comparisons", [])
        ):
            print("Using exported configuration-comparison data to switch the active vessel configuration without rerunning MATLAB.")
            return

    analysis_script = resolve_path(
        config_path.parent,
        str(auto_cfg.get("analysis_script", "..\\analyze_vermicomposter_heater.py")),
    )
    command = [sys.executable, str(analysis_script), "--study-config", str(config_path)]

    if auto_cfg.get("output_dir"):
        command.extend(["--output-dir", str(resolve_path(config_path.parent, str(auto_cfg["output_dir"])))])
    if auto_cfg.get("model_grid_points") is not None:
        command.extend(["--model-grid-points", str(int(auto_cfg["model_grid_points"]))])
    if auto_cfg.get("plot_points") is not None:
        command.extend(["--plot-points", str(int(auto_cfg["plot_points"]))])
    if auto_cfg.get("matlab"):
        command.extend(["--matlab", str(Path(str(auto_cfg["matlab"])).expanduser())])

    print("Detected model_overrides change; refreshing MATLAB export before plotting...")
    subprocess.run(command, cwd=str(analysis_script.parent), check=True)


def apply_config_overrides(payload: dict, config: dict) -> dict:
    payload = json.loads(json.dumps(payload))

    if "recommended_definition" in config:
        payload["meta"]["recommended_definition"] = config["recommended_definition"]
    if "constraint_plot_definition" in config:
        payload["meta"]["constraint_plot_definition"] = config["constraint_plot_definition"]

    summary_inputs = config.get("summary_inputs", {})
    bed_setpoint_c = float(summary_inputs.get("bedSetpoint_C", 20.0))
    current_overrides = normalized_model_overrides(config)
    env_inputs = effective_environment_inputs(summary_inputs, current_overrides)
    greenhouse_air_c = float(env_inputs.get("greenhouseAir_C", -15.0))
    max_air_outlet_c = float(
        current_overrides.get("limits", {}).get(
            "maxAirOutletTemp_C",
            payload.get("limits", {}).get("maxAirOutletTemp_C", 30.0),
        )
    )
    target_duty = configured_target_duty(payload, config)
    if not isinstance(payload.get("sweep"), dict):
        payload["sweep"] = {}
    payload["sweep"]["targetDuty"] = target_duty
    if "requiredHeat_W" in payload:
        payload["minimumFlowAtTargetDuty_Lpm"] = minimum_total_flow_for_outlet_limit_python(
            float(payload["requiredHeat_W"]),
            target_duty,
            bed_setpoint_c,
            greenhouse_air_c,
            max_air_outlet_c,
        )

    overrides = config.get("wire_study_target_overrides", {})
    if overrides:
        override_by_label = {item["label"]: item for item in payload["wire_diameter_sweep"]}
        for label, target in overrides.items():
            if label not in override_by_label:
                continue
            entry = override_by_label[label]
            if "power_W" in target:
                entry["targetPower_W"] = target["power_W"]
            if "flow_Lpm" in target:
                entry["targetFlow_Lpm"] = target["flow_Lpm"]
            if "current_A" in target:
                entry["targetCurrent_A"] = target["current_A"]

    desired_label = requested_active_configuration_label(config)
    for item in payload.get("configuration_comparisons", []):
        if "Qreq_W" in item:
            item["flowTargetDuty_Lpm"] = minimum_total_flow_for_outlet_limit_python(
                float(item["Qreq_W"]),
                target_duty,
                bed_setpoint_c,
                greenhouse_air_c,
                max_air_outlet_c,
            )
    if desired_label:
        for idx, item in enumerate(payload.get("configuration_comparisons", []), start=1):
            if same_configuration_label(desired_label, item.get("label", "")):
                payload["active_configuration"] = {
                    "index": idx,
                    "label": item["label"],
                    "color": item.get("color", [0.15, 0.15, 0.15]),
                }
                payload["design_points"] = item.get("design_points", [])
                payload["requiredHeat_W"] = item.get("Qreq_W", payload.get("requiredHeat_W"))
                payload["minimumFlowAt100Duty_Lpm"] = item.get("flow100_Lpm", payload.get("minimumFlowAt100Duty_Lpm"))
                payload["minimumFlowAtTargetDuty_Lpm"] = item.get("flowTargetDuty_Lpm", payload.get("minimumFlowAtTargetDuty_Lpm"))
                break

    return payload


def requested_active_configuration_label(config: dict) -> str | None:
    bin_overrides = normalized_model_overrides(config).get("bin", {})
    if not isinstance(bin_overrides, dict) or not bin_overrides:
        return None
    open_top = bool(bin_overrides.get("openTop", True))
    vent_count = int(bin_overrides.get("ventHoleCount", 0))
    vent_diameter_m = float(bin_overrides.get("ventHoleDiameter_m", 0.0))
    if open_top:
        return "Uncovered top"
    if vent_count > 0 and vent_diameter_m > 0:
        return "Covered top + vent"
    return "Covered top"


def configured_target_duty(payload: dict, config: dict) -> float:
    sweep_cfg = normalized_model_overrides(config).get("sweep", {})
    if isinstance(sweep_cfg, dict) and "targetDuty" in sweep_cfg:
        return float(sweep_cfg["targetDuty"])
    return float(payload.get("sweep", {}).get("targetDuty", 0.20))


def minimum_total_flow_for_outlet_limit_python(
    q_required_w: float,
    duty_fraction: float,
    bed_setpoint_c: float,
    greenhouse_air_c: float,
    max_air_outlet_c: float,
    pressure_pa: float = 101325.0,
) -> float:
    delta_t_use_k = max_air_outlet_c - bed_setpoint_c
    if delta_t_use_k <= 0 or duty_fraction <= 0:
        return float("inf")
    props = air_props(greenhouse_air_c, pressure_pa)
    mdot_min_kg_s = (q_required_w / duty_fraction) / max(props["cp_J_kgK"] * delta_t_use_k, 1e-12)
    return mdot_min_kg_s / max(props["rho_kg_m3"], 1e-12) * 60.0 * 1000.0


def effective_limits(payload: dict, config: dict) -> dict:
    limits = dict(payload.get("limits", {}))
    cfg_limits = normalized_model_overrides(config).get("limits", {})
    if isinstance(cfg_limits, dict):
        limits.update(cfg_limits)
    limits.setdefault("maxHoleVelocity_m_s", np.inf)
    return limits


def export_consistency_warnings(payload: dict, config: dict) -> list[str]:
    warnings: list[str] = []
    current_overrides = normalized_model_overrides(config)
    sweep_cfg = current_overrides.get("sweep", {})
    payload_sweep = payload.get("sweep", {})
    export_snapshot = payload.get("meta", {}).get("model_overrides_snapshot")

    def compare_range(label: str, cfg_values: list[float] | None, payload_key: str) -> None:
        if not cfg_values or payload_key not in payload_sweep:
            return
        payload_values = payload_sweep.get(payload_key, [])
        if not isinstance(payload_values, list) or not payload_values:
            return
        if abs(float(cfg_values[0]) - float(payload_values[0])) > 1e-9 or abs(float(cfg_values[-1]) - float(payload_values[-1])) > 1e-9 or len(cfg_values) != len(payload_values):
            warnings.append(
                f"Configured {label} sweep ({cfg_values[0]:.1f} to {cfg_values[-1]:.1f}, {len(cfg_values)} points) "
                f"differs from loaded export ({payload_values[0]:.1f} to {payload_values[-1]:.1f}, {len(payload_values)} points); "
                "rerun analyze_vermicomposter_heater.py to update the physics dataset."
            )

    compare_range("voltage", sweep_array_from_config(sweep_cfg, "voltage_V", "voltage_V"), "voltage_V")
    compare_range("flow", sweep_array_from_config(sweep_cfg, "totalFlow_Lpm", "totalFlow_Lpm"), "totalFlow_Lpm")
    compare_range("wire-diameter", sweep_array_from_config(sweep_cfg, "wireDiameter_mm", "wireDiameter_mm"), "wireDiameter_mm")
    compare_range("coil-diameter study", sweep_array_from_config(sweep_cfg, "coilMeanD_mm", "coilMeanD_mm"), "coilMeanD_mm")
    compare_range("wire-length study", sweep_array_from_config(sweep_cfg, "wireLength_m", "wireLength_m"), "wireLength_m")
    if export_snapshot and not can_switch_active_bin_from_export(export_snapshot, current_overrides):
        if export_snapshot != current_overrides:
            warnings.append(
                "Current model_overrides physics inputs differ from the loaded export beyond a pure top-mode switch; "
                "rerun analyze_vermicomposter_heater.py so vessel geometry, heater geometry, and other physics-side changes "
                "are reflected in the dataset."
            )
    return warnings


def wire_series_field(entry: dict, preferred: str, fallback: str) -> list:
    values = entry.get(preferred)
    if isinstance(values, list) and any(value is not None for value in values):
        return values
    return entry[fallback]


def active_reference_point(payload: dict) -> dict:
    rec = payload.get("recommended", {})
    if isinstance(rec, dict) and rec.get("pitch_mm") is not None:
        return rec
    best = payload.get("best_available", {})
    if isinstance(best, dict) and best.get("pitch_mm") is not None:
        return best
    return {}


def build_active_grids(payload: dict) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray]:
    points = payload["design_points"]
    rec = active_reference_point(payload)
    if not rec:
        raise ValueError("No recommended or best-available reference point is available for plotting.")
    pitch = float(rec["pitch_mm"])
    filtered = [pt for pt in points if abs(float(pt["pitch_mm"]) - pitch) < 1e-9]

    voltages = np.array(sorted({float(pt["voltage_V"]) for pt in filtered}), dtype=float)
    flows = np.array(sorted({float(pt["totalFlow_Lpm"]) for pt in filtered}), dtype=float)
    shape = (flows.size, voltages.size)

    grids = {
        "QtoBed_W": np.full(shape, np.nan, dtype=float),
        "totalCurrent_A": np.full(shape, np.nan, dtype=float),
        "airOutlet_C": np.full(shape, np.nan, dtype=float),
        "wireMax_C": np.full(shape, np.nan, dtype=float),
        "deltaP_Pa": np.full(shape, np.nan, dtype=float),
        "holeVelocity_m_s": np.full(shape, np.nan, dtype=float),
        "waterLoss_kg_day": np.full(shape, np.nan, dtype=float),
        "latentEvap_W": np.full(shape, np.nan, dtype=float),
    }

    for pt in filtered:
        row = int(np.where(flows == float(pt["totalFlow_Lpm"]))[0][0])
        col = int(np.where(voltages == float(pt["voltage_V"]))[0][0])
        grids["QtoBed_W"][row, col] = float(pt["QtoBed_W"])
        grids["totalCurrent_A"][row, col] = float(pt["totalCurrent_A"])
        grids["airOutlet_C"][row, col] = float(pt["airOutlet_C"])
        grids["wireMax_C"][row, col] = float(pt["wireMax_C"])
        grids["deltaP_Pa"][row, col] = float(pt["deltaP_Pa"])
        grids["holeVelocity_m_s"][row, col] = float(pt.get("holeVelocity_m_s", np.nan))
        grids["waterLoss_kg_day"][row, col] = float(pt.get("waterLoss_kg_day", np.nan))
        grids["latentEvap_W"][row, col] = float(pt.get("latentEvap_W", np.nan))

    return grids, voltages, flows


def resample_regular_grid(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, nx: int, ny: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x_new = np.linspace(x.min(), x.max(), nx)
    y_new = np.linspace(y.min(), y.max(), ny)
    z_x = np.vstack([np.interp(x_new, x, row) for row in z])
    z_xy = np.vstack([np.interp(y_new, y, z_x[:, j]) for j in range(z_x.shape[1])]).T
    return x_new, y_new, z_xy


def contour_levels(z: np.ndarray, count: int) -> np.ndarray:
    finite = z[np.isfinite(z)]
    if finite.size == 0:
        return np.array([], dtype=float)
    z_min = float(np.min(finite))
    z_max = float(np.max(finite))
    if np.isclose(z_min, z_max):
        return np.array([z_min], dtype=float)
    return np.unique(np.round(np.linspace(z_min, z_max, count), 10))


def plot_constraint_maps(payload: dict, output_dir: Path, grid_points: int, curve_count: int) -> None:
    grids, voltages, flows = build_active_grids(payload)
    rec = active_reference_point(payload)
    point_label = "recommended" if payload.get("recommended") else "best-available reference"

    panels = [
        ("QtoBed_W", "Heat to Bed (W)"),
        ("airOutlet_C", "Outlet Air Temperature (C)"),
        ("wireMax_C", "Max Wire Temperature (C)"),
        ("deltaP_Pa", "Total Pressure Drop (Pa)"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.subplots_adjust(left=0.07, right=0.97, top=0.88, bottom=0.14, hspace=0.26, wspace=0.22)
    fig.suptitle(
        f"Constraint and Performance Maps at {float(rec['pitch_mm']):.1f} mm Pitch\n"
        f"Pitch taken from the active {point_label} point; each panel shows {curve_count} contour curves",
        fontsize=15,
        fontweight="bold",
    )

    for ax, (field, title) in zip(axes.flat, panels):
        x_new, y_new, z_new = resample_regular_grid(voltages, flows, grids[field], grid_points, grid_points)
        x_mesh, y_mesh = np.meshgrid(x_new, y_new)
        levels = contour_levels(z_new, curve_count)
        contour = ax.contour(x_mesh, y_mesh, z_new, levels=levels, cmap="viridis", linewidths=1.5)
        if contour.levels.size > 0:
            ax.clabel(contour, inline=True, fontsize=7, fmt="%.0f")
        ax.set_xlabel("Voltage per branch (V)")
        ax.set_ylabel("Total flow (L/min)")
        ax.set_title(title)
        ax.grid(True, alpha=0.15)

    fig.text(0.02, 0.05, payload["meta"]["constraint_plot_definition"], fontsize=9, wrap=True)
    fig.savefig(output_dir / "constraint_performance_maps.png", dpi=220)
    plt.close(fig)


def plot_total_current_curves(payload: dict, output_dir: Path, curve_count: int) -> None:
    grids, voltages, flows = build_active_grids(payload)
    current_grid = grids["totalCurrent_A"]

    if voltages.size == 0 or flows.size == 0:
        return

    chosen = np.unique(np.linspace(0, flows.size - 1, min(curve_count, flows.size), dtype=int))
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(chosen)))

    fig, ax = plt.subplots(figsize=(10, 6.5))
    fig.subplots_adjust(left=0.10, right=0.97, top=0.88, bottom=0.14)
    fig.suptitle(
        "Total Current (A)\n"
        "Current versus voltage for selected total flow rates",
        fontsize=15,
        fontweight="bold",
    )

    for color, flow_idx in zip(colors, chosen):
        ax.plot(
            voltages,
            current_grid[flow_idx, :],
            color=color,
            linewidth=1.8,
            label=f"{flows[flow_idx]:.1f} L/min",
        )

    ax.set_xlabel("Voltage per branch (V)")
    ax.set_ylabel("Total current (A)")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=8, ncols=2, title="Total flow")
    fig.savefig(output_dir / "total_current_curves.png", dpi=220)
    plt.close(fig)


def plot_moisture_maps(payload: dict, output_dir: Path, grid_points: int, curve_count: int) -> None:
    grids, voltages, flows = build_active_grids(payload)
    rec = active_reference_point(payload)
    point_label = "recommended" if payload.get("recommended") else "best-available reference"

    panels = [
        ("holeVelocity_m_s", "Hole Velocity (m/s)"),
        ("waterLoss_kg_day", "Water Loss (kg/24 h)"),
        ("latentEvap_W", "Latent Evaporation Load (W)"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.6))
    fig.subplots_adjust(left=0.06, right=0.98, top=0.86, bottom=0.18, wspace=0.24)
    fig.suptitle(
        f"Perforation and Moisture Maps at {float(rec['pitch_mm']):.1f} mm Pitch\n"
        f"Pitch taken from the active {point_label} point",
        fontsize=15,
        fontweight="bold",
    )

    evap_limit = float(payload["limits"].get("maxWaterLoss_kg_day", np.nan))
    x_loss, y_loss, z_loss = resample_regular_grid(
        voltages, flows, grids["waterLoss_kg_day"], grid_points, grid_points
    )
    x_loss_mesh, y_loss_mesh = np.meshgrid(x_loss, y_loss)

    for ax, (field, title) in zip(axes.flat, panels):
        x_new, y_new, z_new = resample_regular_grid(voltages, flows, grids[field], grid_points, grid_points)
        x_mesh, y_mesh = np.meshgrid(x_new, y_new)
        levels = contour_levels(z_new, curve_count)
        contour = ax.contour(x_mesh, y_mesh, z_new, levels=levels, cmap="cividis", linewidths=1.5)
        if contour.levels.size > 0:
            ax.clabel(contour, inline=True, fontsize=7, fmt="%.1f")
        if np.isfinite(evap_limit):
            ax.contour(
                x_loss_mesh,
                y_loss_mesh,
                z_loss,
                levels=[evap_limit],
                colors="red",
                linewidths=1.6,
                linestyles="--",
            )
        ax.set_xlabel("Voltage per branch (V)")
        ax.set_ylabel("Total flow (L/min)")
        ax.set_title(title)
        ax.grid(True, alpha=0.15)

    fig.text(
        0.02,
        0.06,
        "Hole velocity is computed from released flow divided by total perforation area. "
        "Water loss is the modeled sprinkler make-up requirement; the dashed red contour is the configured 24 h water-loss "
        "limit boundary, so its position on the hole-velocity panel shows the condition-dependent allowable velocity envelope. "
        "The wetted area in this study is derived from the Frederickson et al. (2007) Table 3 particle-size distribution "
        "unless the JSON switches to another basis.",
        fontsize=9,
        wrap=True,
    )
    fig.savefig(output_dir / "moisture_perforation_maps.png", dpi=220)
    plt.close(fig)


def plot_summer_cooling_curves(payload: dict, config: dict, output_dir: Path) -> None:
    points = sorted(rank_cooling_points(payload, config), key=lambda pt: float(pt["totalFlow_Lpm"]))
    if not points:
        return

    flows = np.array([float(pt["totalFlow_Lpm"]) for pt in points], dtype=float)
    tb = np.array([float(pt["TbottomFullPower_C"]) for pt in points], dtype=float)
    tt = np.array([float(pt["TtopFullPower_C"]) for pt in points], dtype=float)
    tair = np.array([float(pt["airOutlet_C"]) for pt in points], dtype=float)
    tin = np.array([float(pt.get("airInlet_C", np.nan)) for pt in points], dtype=float)
    qcool = np.array([max(-float(pt["QtoBed_W"]), 0.0) for pt in points], dtype=float)
    water = np.array([float(pt.get("waterLoss_kg_day", np.nan)) for pt in points], dtype=float)
    spot_load = np.array([float(pt.get("spotCoolerLoad_W", np.nan)) for pt in points], dtype=float)
    spot_power = np.array([float(pt.get("spotCoolerPower_W", np.nan)) for pt in points], dtype=float)
    blower_power = np.array([float(pt.get("assistBlowerPower_W", np.nan)) for pt in points], dtype=float)
    blower_avail = np.array([float(pt.get("assistBlowerAvailablePressure_Pa", np.nan)) for pt in points], dtype=float)
    spot_cond = np.array([float(pt.get("spotCoolerCondensate_kg_day", np.nan)) for pt in points], dtype=float)
    dp = np.array([float(pt["deltaP_Pa"]) for pt in points], dtype=float)
    hole_u = np.array([float(pt.get("holeVelocity_m_s", np.nan)) for pt in points], dtype=float)

    cooling_cfg = cooling_mode_config(config)
    opt = cooling_cfg["optimization"]
    limits = cooling_cfg["limits"]
    required_cooling = float(payload.get("requiredCooling_W", np.nan))
    spot_cfg = dict(cooling_cfg.get("spot_cooler", {}))

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.subplots_adjust(left=0.08, right=0.96, top=0.90, bottom=0.12, hspace=0.28, wspace=0.28)
    fig.suptitle(
        "Summer Spot-Cooler + Assist-Blower Performance Sweep\n"
        "The spot cooler sets plenum air temperature; the downstream assist blower sets deliverable flow and static pressure",
        fontsize=15,
        fontweight="bold",
    )

    ax = axes[0, 0]
    ax.plot(flows, tb, label="Bottom node", linewidth=2.0)
    ax.plot(flows, tt, label="Top node", linewidth=2.0)
    ax.axhline(float(opt["max_bottom_temp_C"]), color="tab:red", linestyle="--", linewidth=1.2, label="Bottom limit")
    ax.axhline(float(opt["max_top_temp_C"]), color="tab:orange", linestyle=":", linewidth=1.2, label="Top limit")
    ax.set_xlabel("Total flow (L/min)")
    ax.set_ylabel("Bed temperature (C)")
    ax.set_title("Bed Temperatures")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=8)

    ax = axes[0, 1]
    ax.plot(flows, qcool, color="tab:blue", linewidth=2.0, label="Net bed-cooling capacity")
    if np.isfinite(required_cooling):
        ax.axhline(required_cooling, color="tab:red", linestyle="--", linewidth=1.2, label="Lumped cooling requirement")
    ax2 = ax.twinx()
    ax2.plot(flows, tin, color="tab:purple", linewidth=1.4, linestyle="--", label="Spot-cooler supply air")
    ax2.plot(flows, tair, color="tab:green", linewidth=1.8, label="Aeration outlet air")
    ax.set_xlabel("Total flow (L/min)")
    ax.set_ylabel("Cooling capacity (W)")
    ax2.set_ylabel("Outlet air temperature (C)")
    ax.set_title("Cooling Capacity and Outlet Air")
    ax.grid(True, alpha=0.25)
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [ln.get_label() for ln in lines], loc="best", fontsize=8)

    ax = axes[1, 0]
    ax.plot(flows, water, color="tab:cyan", linewidth=2.0, label="Substrate water loss")
    if np.isfinite(spot_cond).any():
        ax.plot(flows, spot_cond, color="tab:blue", linewidth=1.6, linestyle="--", label="Spot-cooler condensate")
    if np.isfinite(float(limits.get("maxWaterLoss_kg_day", np.nan))):
        ax.axhline(float(limits["maxWaterLoss_kg_day"]), color="tab:red", linestyle="--", linewidth=1.2, label="Water-loss limit")
    ax2 = ax.twinx()
    ax2.plot(flows, spot_load, color="tab:purple", linewidth=1.8, label="Spot-cooler sensible load")
    if np.isfinite(spot_power).any():
        ax2.plot(flows, spot_power, color="tab:pink", linewidth=1.4, linestyle=":", label="Spot-cooler power")
    if np.isfinite(blower_power).any():
        ax2.plot(flows, blower_power, color="tab:gray", linewidth=1.4, linestyle="-.", label="Assist-blower power")
    ax.set_xlabel("Total flow (L/min)")
    ax.set_ylabel("Water loss (kg/24 h)")
    ax2.set_ylabel("Cooling-system load / power (W)")
    ax.set_title("Water Loss and Cooling-System Power")
    ax.grid(True, alpha=0.25)
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [ln.get_label() for ln in lines], loc="best", fontsize=8)

    ax = axes[1, 1]
    ax.plot(flows, dp, color="tab:brown", linewidth=2.0, label="Pressure drop")
    if np.isfinite(blower_avail).any():
        ax.plot(flows, blower_avail, color="tab:blue", linewidth=1.6, linestyle="--", label="Assist-blower available static")
    if np.isfinite(float(limits.get("maxPressureDrop_Pa", np.nan))):
        ax.axhline(float(limits["maxPressureDrop_Pa"]), color="tab:red", linestyle="--", linewidth=1.2, label="dP limit")
    ax2 = ax.twinx()
    ax2.plot(flows, hole_u, color="tab:olive", linewidth=1.8, label="Hole velocity")
    if np.isfinite(float(limits.get("maxHoleVelocity_m_s", np.inf))):
        ax2.axhline(float(limits["maxHoleVelocity_m_s"]), color="tab:orange", linestyle=":", linewidth=1.2, label="Velocity limit")
    ax.set_xlabel("Total flow (L/min)")
    ax.set_ylabel("Pressure drop (Pa)")
    ax2.set_ylabel("Hole velocity (m/s)")
    ax.set_title("Pressure Drop and Hole Velocity")
    ax.grid(True, alpha=0.25)
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [ln.get_label() for ln in lines], loc="best", fontsize=8)

    fig.text(
        0.02,
        0.03,
        "The summer sweep holds electrical heater power at zero. A Whynter spot cooler first conditions an insulated plenum. "
        "A downstream assist blower then draws from that plenum and pushes the cooled air through the aeration network. "
        "The aeration tube applies wall exchange, direct jet exchange, and substrate-evaporation losses segment-by-segment. "
        "Negative Q_to_bed denotes net cooling of the substrate.",
        fontsize=9,
        wrap=True,
    )
    fig.savefig(output_dir / "summer_spot_cooler_curves.png", dpi=220)
    plt.close(fig)


def compute_radial_profile(payload: dict, config: dict) -> dict:
    radial_cfg = dict(config.get("radial_profile", {}))
    if not radial_cfg.get("enabled", True):
        return {}

    rec = active_reference_point(payload)
    if not rec:
        return {}

    model_overrides = normalized_model_overrides(config)
    summary_inputs = config.get("summary_inputs", {})
    mode = operating_mode(config)
    env_inputs = effective_environment_inputs(summary_inputs, model_overrides)
    if mode == COOLING_MODE:
        cooling_cfg = cooling_mode_config(config)
        env_inputs["greenhouseAir_C"] = float(cooling_cfg["ambientAir_C"])
        env_inputs["pressure_Pa"] = float(cooling_cfg.get("pressure_Pa", env_inputs.get("pressure_Pa", 101325.0)))
    bin_inputs = effective_bin_inputs(summary_inputs, model_overrides)
    aer_inputs = effective_aeration_inputs(summary_inputs, model_overrides)
    pressure_Pa = float(env_inputs.get("pressure_Pa", 101325.0))
    greenhouse_air_C = float(env_inputs.get("greenhouseAir_C", -15.0))
    rec_air_C = float(rec.get("airOutlet_C", np.nan))
    if not np.isfinite(rec_air_C):
        return {}

    ref_mode = str(radial_cfg.get("bed_reference", "bottom_node")).strip().lower()
    T_bottom = float(rec.get("TbottomFullPower_C", np.nan))
    T_top = float(rec.get("TtopFullPower_C", np.nan))
    if ref_mode == "mean_nodes" and np.isfinite(T_bottom) and np.isfinite(T_top):
        bed_ref_C = 0.5 * (T_bottom + T_top)
    elif ref_mode == "colder_node" and np.isfinite(T_bottom) and np.isfinite(T_top):
        bed_ref_C = min(T_bottom, T_top)
    elif np.isfinite(T_bottom):
        bed_ref_C = T_bottom
    elif np.isfinite(T_top):
        bed_ref_C = T_top
    else:
        bed_ref_C = float(summary_inputs.get("bedSetpoint_C", 20.0))

    evaporation_inputs = model_overrides.get("evaporation", {})
    surface_set_C = float(evaporation_inputs.get("surfaceTemp_C", summary_inputs.get("bedSetpoint_C", bed_ref_C)))
    sink_C = min(max(surface_set_C, bed_ref_C), rec_air_C)
    wetted = derived_wetted_area_details(summary_inputs, model_overrides)
    fill_volume_m3 = float(wetted.get("fill_volume_m3", np.nan))
    wetted_area_m2 = float(wetted.get("area_m2", np.nan))
    if not (np.isfinite(fill_volume_m3) and fill_volume_m3 > 0 and np.isfinite(wetted_area_m2) and wetted_area_m2 > 0):
        return {}

    area_density_m2_m3 = wetted_area_m2 / fill_volume_m3
    aeration_overrides = model_overrides.get("aeration", {})
    calibration = model_overrides.get("calibration", {})
    h_sensible = float(aeration_overrides.get("bedH_W_m2K", 0.0)) * float(calibration.get("bedHTMultiplier", 1.0))
    if h_sensible <= 0:
        return {}

    n_tubes = max(1.0, float(aer_inputs.get("nParallelTubes", 1.0)))
    aer_len_m = max(1e-9, float(aer_inputs.get("length_m", 1.0)))
    release_fraction = float(aeration_overrides.get("releaseFraction", 1.0))
    q_total_m3_s = float(rec.get("totalFlow_Lpm", 0.0)) / 60000.0
    air_in = air_props(greenhouse_air_C, pressure_Pa)
    mdot_tube_kg_s = air_in["rho_kg_m3"] * q_total_m3_s / n_tubes
    mdot_release_per_length_kg_s_m = release_fraction * mdot_tube_kg_s / aer_len_m
    if mdot_release_per_length_kg_s_m <= 0:
        return {}

    aer_od_m = float(aer_inputs.get("OD_mm", 44.0)) / 1000.0
    r0_m = 0.5 * aer_od_m
    max_distance_m = float(radial_cfg.get("max_distance_m", 0.20))
    point_count = max(25, int(radial_cfg.get("point_count", 300)))
    threshold_C = float(radial_cfg.get("temperature_threshold_C", 27.0))
    spread_rate_m_per_m = max(0.0, float(radial_cfg.get("spread_rate_m_per_m", 0.15)))

    initial_mode = str(radial_cfg.get("initial_plume_radius_mode", "hole_radius")).strip().lower()
    initial_override = radial_cfg.get("initial_plume_radius_m")
    hole_diameter_m = float(model_overrides.get("perforation", {}).get("holeDiameter_m", 0.005))
    if initial_override is not None:
        plume_radius0_m = max(1e-9, float(initial_override))
    elif initial_mode == "hole_diameter":
        plume_radius0_m = max(1e-9, hole_diameter_m)
    else:
        plume_radius0_m = max(1e-9, 0.5 * hole_diameter_m)

    max_plume_radius_raw = radial_cfg.get("max_plume_radius_m")
    max_plume_radius_m = np.inf if max_plume_radius_raw is None else max(plume_radius0_m, float(max_plume_radius_raw))
    distances_m = np.linspace(0.0, max_distance_m, point_count)
    plume_radius_m = np.minimum(max_plume_radius_m, plume_radius0_m + spread_rate_m_per_m * distances_m)
    temperatures_C = np.full_like(distances_m, rec_air_C, dtype=float)
    bulk_vapor_density_kg_m3 = np.full_like(distances_m, np.nan, dtype=float)
    local_relative_humidity = np.full_like(distances_m, np.nan, dtype=float)
    sensible_load_W_m = np.zeros_like(distances_m, dtype=float)
    latent_load_W_m = np.zeros_like(distances_m, dtype=float)
    local_h_mass_m_s = np.zeros_like(distances_m, dtype=float)
    local_h_latent_W_m2K = np.zeros_like(distances_m, dtype=float)
    contact_area_per_length_m2_m = np.pi * area_density_m2_m3 * plume_radius_m**2

    evaporation_inputs = model_overrides.get("evaporation", {})
    relative_humidity = float(evaporation_inputs.get("relativeHumidity", 0.8))
    lewis_factor = float(evaporation_inputs.get("lewisFactor", 1.0))
    evaporation_multiplier = float(calibration.get("evaporationHMultiplier", 1.0))
    rho_v_inlet = relative_humidity * water_vapor_density_at_saturation_kg_m3(greenhouse_air_C + 273.15)
    rho_v_bulk = min(rho_v_inlet, water_vapor_density_at_saturation_kg_m3(rec_air_C + 273.15))
    bulk_vapor_density_kg_m3[0] = rho_v_bulk
    local_relative_humidity[0] = min(
        1.0,
        rho_v_bulk / max(water_vapor_density_at_saturation_kg_m3(rec_air_C + 273.15), 1e-12),
    )

    for idx in range(distances_m.size - 1):
        dx_m = float(distances_m[idx + 1] - distances_m[idx])
        if dx_m <= 0:
            continue

        T_local_C = float(temperatures_C[idx])
        if T_local_C <= sink_C + 1e-12:
            temperatures_C[idx + 1] = sink_C
            bulk_vapor_density_kg_m3[idx + 1] = rho_v_bulk
            local_relative_humidity[idx + 1] = min(
                1.0,
                rho_v_bulk / max(water_vapor_density_at_saturation_kg_m3(sink_C + 273.15), 1e-12),
            )
            continue

        props = air_props(T_local_C, pressure_Pa)
        local_surface_C = min(max(surface_set_C, bed_ref_C), T_local_C)
        local_surface_K = local_surface_C + 273.15
        rho_v_surface = water_vapor_density_at_saturation_kg_m3(local_surface_K)
        driving_kg_m3 = max(rho_v_surface - rho_v_bulk, 0.0)

        h_mass_m_s = evaporation_multiplier * h_sensible / max(
            props["rho_kg_m3"] * props["cp_J_kgK"] * lewis_factor ** (2.0 / 3.0),
            1e-12,
        )
        local_h_mass_m_s[idx] = h_mass_m_s

        area_contact = float(contact_area_per_length_m2_m[idx])
        q_sensible_prime_W_m = h_sensible * area_contact * max(T_local_C - sink_C, 0.0)
        m_evap_prime_kg_s_m = h_mass_m_s * area_contact * driving_kg_m3

        vdot_release_per_length_m2_s = mdot_release_per_length_kg_s_m / max(props["rho_kg_m3"], 1e-12)
        rho_v_cap_kg_m3 = water_vapor_density_at_saturation_kg_m3(T_local_C + 273.15)
        m_capacity_step_kg_s = max(vdot_release_per_length_m2_s * (rho_v_cap_kg_m3 - rho_v_bulk), 0.0)
        m_evap_step_kg_s = min(m_evap_prime_kg_s_m * dx_m, m_capacity_step_kg_s)

        hfg_J_kg = latent_heat_vaporization_water_J_kg(local_surface_K)
        q_sensible_step_W = q_sensible_prime_W_m * dx_m
        q_latent_step_W = m_evap_step_kg_s * hfg_J_kg
        available_step_W = mdot_release_per_length_kg_s_m * props["cp_J_kgK"] * max(T_local_C - sink_C, 0.0)
        total_step_W = q_sensible_step_W + q_latent_step_W
        if total_step_W > available_step_W and total_step_W > 0:
            scale = available_step_W / total_step_W
            q_sensible_step_W *= scale
            q_latent_step_W *= scale
            m_evap_step_kg_s *= scale
        else:
            scale = 1.0

        sensible_load_W_m[idx] = q_sensible_step_W / max(dx_m, 1e-12)
        latent_load_W_m[idx] = q_latent_step_W / max(dx_m, 1e-12)
        local_h_latent_W_m2K[idx] = latent_load_W_m[idx] / max(area_contact * max(T_local_C - sink_C, 1e-9), 1e-12)

        T_next_C = T_local_C - (q_sensible_step_W + q_latent_step_W) / max(
            mdot_release_per_length_kg_s_m * props["cp_J_kgK"],
            1e-12,
        )
        T_next_C = max(T_next_C, sink_C)
        temperatures_C[idx + 1] = T_next_C

        rho_v_bulk = rho_v_bulk + m_evap_step_kg_s / max(vdot_release_per_length_m2_s, 1e-12)
        rho_v_bulk = min(rho_v_bulk, water_vapor_density_at_saturation_kg_m3(T_next_C + 273.15))
        bulk_vapor_density_kg_m3[idx + 1] = rho_v_bulk
        local_relative_humidity[idx + 1] = min(
            1.0,
            rho_v_bulk / max(water_vapor_density_at_saturation_kg_m3(T_next_C + 273.15), 1e-12),
        )

    if distances_m.size > 1:
        sensible_load_W_m[-1] = sensible_load_W_m[-2]
        latent_load_W_m[-1] = latent_load_W_m[-2]
        local_h_mass_m_s[-1] = local_h_mass_m_s[-2]
        local_h_latent_W_m2K[-1] = local_h_latent_W_m2K[-2]

    h_latent = float(np.nanmean(local_h_latent_W_m2K)) if np.isfinite(np.nanmean(local_h_latent_W_m2K)) else 0.0
    h_effective = h_sensible + h_latent
    sink_strength = float(np.nanmean(
        np.pi * area_density_m2_m3 * (h_sensible + np.maximum(local_h_latent_W_m2K, 0.0))
        / max(mdot_release_per_length_kg_s_m * air_in["cp_J_kgK"], 1e-12)
    ))

    safe_distance_m = np.nan
    threshold_state = "not_needed"
    if rec_air_C > threshold_C:
        if sink_C >= threshold_C:
            threshold_state = "not_reached_asymptote"
        else:
            crossing = np.where(temperatures_C <= threshold_C)[0]
            if crossing.size > 0:
                idx = int(crossing[0])
                if idx == 0:
                    safe_distance_m = distances_m[0]
                else:
                    x0 = distances_m[idx - 1]
                    x1 = distances_m[idx]
                    y0 = temperatures_C[idx - 1]
                    y1 = temperatures_C[idx]
                    if np.isclose(y1, y0):
                        safe_distance_m = x1
                    else:
                        frac = (threshold_C - y0) / (y1 - y0)
                        safe_distance_m = x0 + frac * (x1 - x0)
                threshold_state = "reached"
            elif temperatures_C[-1] > threshold_C:
                threshold_state = "beyond_plot_limit"
            else:
                threshold_state = "not_reached_asymptote"

    return {
        "temperatures_C": temperatures_C,
        "distances_m": distances_m,
        "threshold_C": threshold_C,
        "safe_distance_m": safe_distance_m,
        "threshold_state": threshold_state,
        "sink_C": sink_C,
        "bed_ref_C": bed_ref_C,
        "area_density_m2_m3": area_density_m2_m3,
        "h_sensible_W_m2K": h_sensible,
        "h_latent_W_m2K": h_latent,
        "h_effective_W_m2K": h_effective,
        "bulk_vapor_density_kg_m3": bulk_vapor_density_kg_m3,
        "local_relative_humidity": local_relative_humidity,
        "local_h_mass_m_s": local_h_mass_m_s,
        "local_h_latent_W_m2K": local_h_latent_W_m2K,
        "sensible_load_W_m": sensible_load_W_m,
        "latent_load_W_m": latent_load_W_m,
        "total_radial_sensible_W": float(np.trapezoid(sensible_load_W_m, distances_m)),
        "total_radial_latent_W": float(np.trapezoid(latent_load_W_m, distances_m)),
        "initial_bulk_vapor_density_kg_m3": float(bulk_vapor_density_kg_m3[0]),
        "final_bulk_vapor_density_kg_m3": float(bulk_vapor_density_kg_m3[-1]),
        "relative_humidity_inlet": float(local_relative_humidity[0]),
        "relative_humidity_final": float(local_relative_humidity[-1]),
        "evaporation_resolved_locally": True,
        "mdot_release_per_length_kg_s_m": mdot_release_per_length_kg_s_m,
        "plume_radius_m": plume_radius_m,
        "initial_plume_radius_m": plume_radius0_m,
        "spread_rate_m_per_m": spread_rate_m_per_m,
        "max_plume_radius_m": max_plume_radius_m,
        "sink_strength": sink_strength,
        "tube_outer_radius_m": r0_m,
        "point_label": "recommended" if payload.get("recommended") else "best-available reference",
    }


def rule_based_tube_centers(width_m: float, fill_height_m: float, n_tubes: int, layout_cfg: dict) -> tuple[np.ndarray, np.ndarray, str]:
    max_supported = int(layout_cfg.get("maxSupportedTubes", 12))
    if n_tubes < 1 or n_tubes > max_supported:
        raise ValueError(f"rule-based tube layout supports 1 to {max_supported} tubes, got {n_tubes}")

    x_margin = float(layout_cfg.get("wallOffsetWidthFraction", 0.10)) * width_m
    y_margin = float(layout_cfg.get("wallOffsetHeightFraction", 0.10)) * fill_height_m
    x_left = x_margin
    x_right = width_m - x_margin
    x_center = 0.5 * width_m
    x_bottom_row = np.linspace(x_left, x_right, 3)
    y_bottom = y_margin
    y_top = fill_height_m - y_margin
    y_mid = 0.5 * fill_height_m

    def bottom_interior_positions(count: int) -> list[float]:
        if count <= 0:
            return []
        return [float(x) for x in np.linspace(x_left, x_right, count + 2)[1:-1]]

    def side_column_positions(count: int) -> list[float]:
        return [float(y) for y in np.linspace(y_bottom, y_top, count)]

    layouts: dict[int, list[tuple[float, float]]] = {
        1: [(x_center, y_bottom)],
        2: [(x_left, y_bottom), (x_right, y_bottom)],
        3: [(float(x), y_bottom) for x in x_bottom_row],
        4: [(x_left, y_bottom), (x_right, y_bottom), (x_left, y_top), (x_right, y_top)],
        5: [(float(x), y_bottom) for x in x_bottom_row] + [(x_left, y_top), (x_right, y_top)],
        6: [(x_left, float(y)) for y in side_column_positions(3)] + [(x_right, float(y)) for y in side_column_positions(3)],
        7: [
            (x_left, y_bottom),
            (x_right, y_bottom),
            (x_left, y_top),
            (x_right, y_top),
            (x_left, y_mid),
            (x_right, y_mid),
            (x_center, y_bottom),
        ],
        8: [(x_left, y) for y in side_column_positions(3)]
        + [(x_right, y) for y in side_column_positions(3)]
        + [(x, y_bottom) for x in bottom_interior_positions(2)],
        9: [(x_left, y) for y in side_column_positions(3)]
        + [(x_right, y) for y in side_column_positions(3)]
        + [(x, y_bottom) for x in bottom_interior_positions(3)],
        10: [(x_left, y) for y in side_column_positions(3)]
        + [(x_right, y) for y in side_column_positions(3)]
        + [(x, y_bottom) for x in bottom_interior_positions(4)],
        11: [(x_left, y) for y in side_column_positions(4)]
        + [(x_right, y) for y in side_column_positions(4)]
        + [(x, y_bottom) for x in bottom_interior_positions(3)],
        12: [(x_left, y) for y in side_column_positions(4)]
        + [(x_right, y) for y in side_column_positions(4)]
        + [(x, y_bottom) for x in bottom_interior_positions(4)],
    }
    coords = layouts[n_tubes]
    return (
        np.array([xy[0] for xy in coords], dtype=float),
        np.array([xy[1] for xy in coords], dtype=float),
        "rule_based_1to12",
    )


def configured_single_row_tube_center_height_m(layout_cfg: dict, fill_height_m: float, tube_radius_m: float) -> float:
    raw = layout_cfg.get("tubeCenterHeight_m", None)
    if raw not in (None, ""):
        return float(raw)
    return max(float(layout_cfg.get("wallOffsetHeightFraction", 0.10)) * fill_height_m, tube_radius_m)


def compute_tube_layout_and_habitat(payload: dict, config: dict) -> dict:
    layout_cfg = dict(config.get("tube_layout", {}))
    if not layout_cfg.get("enabled", True):
        return {}

    summary_inputs = config.get("summary_inputs", {})
    model_overrides = normalized_model_overrides(config)
    bin_inputs = effective_bin_inputs(summary_inputs, model_overrides)
    aer_inputs = effective_aeration_inputs(summary_inputs, model_overrides)
    if not {"width_m", "height_m", "fillFraction"} <= bin_inputs.keys():
        return {}
    if not {"nParallelTubes", "OD_mm", "length_m"} <= aer_inputs.keys():
        return {}

    width_m = float(bin_inputs["width_m"])
    fill_height_m = float(bin_inputs["height_m"]) * float(bin_inputs["fillFraction"])
    length_m = float(aer_inputs["length_m"])
    n_tubes = int(aer_inputs["nParallelTubes"])
    tube_radius_m = 0.5 * float(aer_inputs["OD_mm"]) / 1000.0

    explicit_centers = layout_cfg.get("tubeCenters_m", [])
    if isinstance(explicit_centers, list) and len(explicit_centers) == n_tubes:
        coords: list[tuple[float, float]] = []
        for item in explicit_centers:
            if isinstance(item, dict) and {"x_m", "y_m"} <= item.keys():
                coords.append((float(item["x_m"]), float(item["y_m"])))
            elif isinstance(item, list) and len(item) == 2:
                coords.append((float(item[0]), float(item[1])))
        if len(coords) == n_tubes:
            x_centers = np.array([xy[0] for xy in coords], dtype=float)
            y_centers = np.array([xy[1] for xy in coords], dtype=float)
            center_mode = "explicit_xy"
        else:
            x_centers = np.array([], dtype=float)
            y_centers = np.array([], dtype=float)
            center_mode = "invalid_explicit_xy"
    else:
        x_centers_raw = layout_cfg.get("tubeCenterX_m", [])
        if isinstance(x_centers_raw, list) and len(x_centers_raw) == n_tubes:
            x_centers = np.array([float(x) for x in x_centers_raw], dtype=float)
            y_center_m = configured_single_row_tube_center_height_m(layout_cfg, fill_height_m, tube_radius_m)
            y_centers = np.full_like(x_centers, y_center_m)
            center_mode = "user-specified-x"
        else:
            layout_mode = str(layout_cfg.get("layoutMode", "rule_based_1to12")).strip().lower()
            if layout_mode in {"rule_based_1to8", "rule_based_1to12"}:
                x_centers, y_centers, center_mode = rule_based_tube_centers(width_m, fill_height_m, n_tubes, layout_cfg)
            else:
                pitch_m = width_m / max(n_tubes, 1)
                x_centers = np.linspace(0.5 * pitch_m, width_m - 0.5 * pitch_m, n_tubes)
                y_center_m = configured_single_row_tube_center_height_m(layout_cfg, fill_height_m, tube_radius_m)
                y_centers = np.full_like(x_centers, y_center_m)
                center_mode = "uniform_across_width"

    radial = compute_radial_profile(payload, config)
    safe_distance_m = float(radial.get("safe_distance_m", np.nan))
    if radial.get("threshold_state") == "not_needed":
        safe_distance_m = 0.0
    exclusion_radius_m = tube_radius_m + max(safe_distance_m, 0.0) if np.isfinite(safe_distance_m) else np.nan

    grid_points = max(250, int(layout_cfg.get("cross_section_grid_points", 700)))
    x_vals = np.linspace(0.0, width_m, grid_points, endpoint=False) + width_m / (2.0 * grid_points)
    y_vals = np.linspace(0.0, fill_height_m, grid_points, endpoint=False) + fill_height_m / (2.0 * grid_points)
    xx, yy = np.meshgrid(x_vals, y_vals)

    pipe_masks = []
    exclusion_masks = []
    for xc, yc in zip(x_centers, y_centers):
        r2 = (xx - xc) ** 2 + (yy - yc) ** 2
        pipe_masks.append(r2 <= tube_radius_m**2)
        if np.isfinite(exclusion_radius_m):
            exclusion_masks.append(r2 <= exclusion_radius_m**2)

    pipe_union = np.logical_or.reduce(pipe_masks) if pipe_masks else np.zeros_like(xx, dtype=bool)
    exclusion_union = np.logical_or.reduce(exclusion_masks) if exclusion_masks else np.zeros_like(xx, dtype=bool)

    wall_masks: list[np.ndarray] = []
    wall_return = {
        "enabled": bool(layout_cfg.get("includeWallConduction", True)),
        "left_depth_m": 0.0,
        "right_depth_m": 0.0,
        "bottom_depth_m": 0.0,
        "left_contactTemp_C": np.nan,
        "right_contactTemp_C": np.nan,
        "bottom_contactTemp_C": np.nan,
        "left_wallTemp_C": np.nan,
        "right_wallTemp_C": np.nan,
        "bottom_wallTemp_C": np.nan,
    }
    wall_inputs = effective_bin_wall_inputs(summary_inputs, model_overrides)
    h_ext = float(summary_inputs.get("externalNaturalConvection_h_W_m2K", np.nan))
    if wall_return["enabled"] and radial and {"sheetThickness_mm", "sheetK_W_mK", "insulationThickness_mm", "insulationK_W_mK"} <= wall_inputs.keys() and np.isfinite(h_ext):
        R_out_m2K_W = (
            float(wall_inputs["sheetThickness_mm"]) / 1000.0 / max(float(wall_inputs["sheetK_W_mK"]), 1e-12)
            + float(wall_inputs["insulationThickness_mm"]) / 1000.0 / max(float(wall_inputs["insulationK_W_mK"]), 1e-12)
            + 1.0 / max(h_ext, 1e-12)
        )
        h_in = max(float(radial.get("h_sensible_W_m2K", 0.0)), 1e-12)
        T_inf = float(summary_inputs.get("greenhouseAir_C", -15.0))
        T_bed = float(radial.get("bed_ref_C", summary_inputs.get("bedSetpoint_C", 20.0)))
        T_safe = float(radial.get("threshold_C", 27.0))
        profile_x = as_float_array(radial.get("distances_m", []))
        profile_T = as_float_array(radial.get("temperatures_C", []))

        def wall_band_depth(wall_temp_C: float, depth_limit_m: float) -> float:
            if not (np.isfinite(wall_temp_C) and depth_limit_m > 0 and wall_temp_C > max(T_bed, T_safe)):
                return 0.0
            if T_bed >= T_safe:
                return depth_limit_m
            return float(np.clip(depth_limit_m * (wall_temp_C - T_safe) / max(wall_temp_C - T_bed, 1e-12), 0.0, depth_limit_m))

        def wall_contact_temperature(distances_to_wall_m: np.ndarray) -> float:
            if distances_to_wall_m.size == 0:
                return float("nan")
            samples = [
                interpolate_profile_value(profile_x, profile_T, max(float(distance), 0.0))
                for distance in distances_to_wall_m
            ]
            finite = [value for value in samples if np.isfinite(value)]
            return float(np.mean(finite)) if finite else float("nan")

        left_mask = np.isclose(x_centers, np.min(x_centers))
        right_mask = np.isclose(x_centers, np.max(x_centers))
        bottom_mask = np.isclose(y_centers, np.min(y_centers))

        left_gap = np.maximum(x_centers[left_mask] - tube_radius_m, 0.0)
        right_gap = np.maximum(width_m - x_centers[right_mask] - tube_radius_m, 0.0)
        bottom_gap = np.maximum(y_centers[bottom_mask] - tube_radius_m, 0.0)

        wall_return["left_contactTemp_C"] = wall_contact_temperature(left_gap)
        wall_return["right_contactTemp_C"] = wall_contact_temperature(right_gap)
        wall_return["bottom_contactTemp_C"] = wall_contact_temperature(bottom_gap)

        for side in ("left", "right", "bottom"):
            contact_temp = wall_return[f"{side}_contactTemp_C"]
            if np.isfinite(contact_temp):
                qpp = (contact_temp - T_inf) / max(1.0 / h_in + R_out_m2K_W, 1e-12)
                wall_return[f"{side}_wallTemp_C"] = contact_temp - qpp / h_in

        left_depth_limit = float(np.min(x_centers[left_mask])) if np.any(left_mask) else 0.0
        right_depth_limit = float(width_m - np.max(x_centers[right_mask])) if np.any(right_mask) else 0.0
        bottom_depth_limit = float(np.min(y_centers[bottom_mask])) if np.any(bottom_mask) else 0.0
        wall_return["left_depth_m"] = wall_band_depth(float(wall_return["left_wallTemp_C"]), left_depth_limit)
        wall_return["right_depth_m"] = wall_band_depth(float(wall_return["right_wallTemp_C"]), right_depth_limit)
        wall_return["bottom_depth_m"] = wall_band_depth(float(wall_return["bottom_wallTemp_C"]), bottom_depth_limit)

        if wall_return["left_depth_m"] > 0:
            wall_masks.append(xx <= wall_return["left_depth_m"])
        if wall_return["right_depth_m"] > 0:
            wall_masks.append(xx >= width_m - wall_return["right_depth_m"])
        if wall_return["bottom_depth_m"] > 0:
            wall_masks.append(yy <= wall_return["bottom_depth_m"])

    wall_union = np.logical_or.reduce(wall_masks) if wall_masks else np.zeros_like(xx, dtype=bool)
    unsafe_union = np.logical_or(exclusion_union, wall_union)

    cross_section_area_m2 = width_m * fill_height_m
    pipe_area_union_m2 = float(np.mean(pipe_union)) * cross_section_area_m2
    exclusion_area_union_m2 = float(np.mean(exclusion_union)) * cross_section_area_m2
    wall_area_union_m2 = float(np.mean(wall_union)) * cross_section_area_m2
    unsafe_area_union_m2 = float(np.mean(unsafe_union)) * cross_section_area_m2
    pipe_volume_m3 = pipe_area_union_m2 * length_m
    exclusion_volume_m3 = exclusion_area_union_m2 * length_m
    wall_volume_m3 = wall_area_union_m2 * length_m
    unsafe_volume_m3 = unsafe_area_union_m2 * length_m
    unsafe_shell_volume_m3 = max(unsafe_volume_m3 - pipe_volume_m3, 0.0)
    fill_volume_m3 = width_m * fill_height_m * length_m
    habitable_volume_m3 = max(fill_volume_m3 - unsafe_volume_m3, 0.0)

    pipe_sum_area_m2 = sum(float(np.mean(mask)) * cross_section_area_m2 for mask in pipe_masks)
    exclusion_sum_area_m2 = sum(float(np.mean(mask)) * cross_section_area_m2 for mask in exclusion_masks)
    pipe_overlap_area_m2 = max(pipe_sum_area_m2 - pipe_area_union_m2, 0.0)
    exclusion_overlap_area_m2 = max(exclusion_sum_area_m2 - exclusion_area_union_m2, 0.0)

    unique_x = np.unique(np.round(x_centers, 12))
    pitches_m = np.diff(unique_x)
    clear_gaps_m = np.maximum(pitches_m - 2.0 * tube_radius_m, 0.0)
    if n_tubes > 1:
        pairwise = []
        for i in range(n_tubes):
            for j in range(i + 1, n_tubes):
                pairwise.append(float(np.hypot(x_centers[i] - x_centers[j], y_centers[i] - y_centers[j])))
        nearest_center_spacing_m = min(pairwise) if pairwise else float("nan")
    else:
        nearest_center_spacing_m = float("nan")
    side_clearance_left_m = float(np.min(x_centers) - tube_radius_m)
    side_clearance_right_m = float(width_m - np.max(x_centers) - tube_radius_m)

    return {
        "width_m": width_m,
        "fill_height_m": fill_height_m,
        "length_m": length_m,
        "n_tubes": n_tubes,
        "tube_radius_m": tube_radius_m,
        "tube_diameter_m": 2.0 * tube_radius_m,
        "center_mode": center_mode,
        "x_centers_m": x_centers,
        "y_centers_m": y_centers,
        "pitches_m": pitches_m,
        "clear_gaps_m": clear_gaps_m,
        "nearest_center_spacing_m": nearest_center_spacing_m,
        "nearest_clear_gap_m": max(nearest_center_spacing_m - 2.0 * tube_radius_m, 0.0) if np.isfinite(nearest_center_spacing_m) else np.nan,
        "side_clearance_left_m": side_clearance_left_m,
        "side_clearance_right_m": side_clearance_right_m,
        "safe_distance_m": safe_distance_m,
        "exclusion_radius_m": exclusion_radius_m,
        "pipe_volume_m3": pipe_volume_m3,
        "exclusion_volume_m3": exclusion_volume_m3,
        "wall_return_volume_m3": wall_volume_m3,
        "unsafe_volume_m3": unsafe_volume_m3,
        "unsafe_shell_volume_m3": unsafe_shell_volume_m3,
        "habitable_volume_m3": habitable_volume_m3,
        "fill_volume_m3": fill_volume_m3,
        "pipe_overlap_volume_m3": pipe_overlap_area_m2 * length_m,
        "exclusion_overlap_volume_m3": exclusion_overlap_area_m2 * length_m,
        "pipe_area_union_m2": pipe_area_union_m2,
        "exclusion_area_union_m2": exclusion_area_union_m2,
        "wall_area_union_m2": wall_area_union_m2,
        "unsafe_area_union_m2": unsafe_area_union_m2,
        "cross_section_grid_points": grid_points,
        "wall_return": wall_return,
    }


def plot_habitat_exclusion_cross_section(payload: dict, config: dict, output_dir: Path) -> None:
    layout = compute_tube_layout_and_habitat(payload, config)
    if not layout:
        return

    fig, ax = plt.subplots(figsize=(9, 6.2))
    fig.subplots_adjust(left=0.08, right=0.97, top=0.88, bottom=0.19)
    ax.set_xlim(0, 1000.0 * layout["width_m"])
    ax.set_ylim(0, 1000.0 * layout["fill_height_m"])
    rect = plt.Rectangle((0, 0), 1000.0 * layout["width_m"], 1000.0 * layout["fill_height_m"],
                         facecolor="#ecfeff", edgecolor="#164e63", linewidth=1.5)
    ax.add_patch(rect)
    wall_return = dict(layout.get("wall_return", {}))
    if float(wall_return.get("left_depth_m", 0.0)) > 0:
        ax.add_patch(
            plt.Rectangle(
                (0, 0),
                1000.0 * float(wall_return["left_depth_m"]),
                1000.0 * layout["fill_height_m"],
                facecolor="#fdba74",
                edgecolor="#ea580c",
                alpha=0.28,
                linewidth=1.0,
            )
        )
    if float(wall_return.get("right_depth_m", 0.0)) > 0:
        ax.add_patch(
            plt.Rectangle(
                (1000.0 * (layout["width_m"] - float(wall_return["right_depth_m"])), 0),
                1000.0 * float(wall_return["right_depth_m"]),
                1000.0 * layout["fill_height_m"],
                facecolor="#fdba74",
                edgecolor="#ea580c",
                alpha=0.28,
                linewidth=1.0,
            )
        )
    if float(wall_return.get("bottom_depth_m", 0.0)) > 0:
        ax.add_patch(
            plt.Rectangle(
                (0, 0),
                1000.0 * layout["width_m"],
                1000.0 * float(wall_return["bottom_depth_m"]),
                facecolor="#fdba74",
                edgecolor="#ea580c",
                alpha=0.28,
                linewidth=1.0,
            )
        )

    for xc, yc in zip(layout["x_centers_m"], layout["y_centers_m"]):
        if np.isfinite(layout["exclusion_radius_m"]):
            excl = plt.Circle((1000.0 * xc, 1000.0 * yc), 1000.0 * layout["exclusion_radius_m"],
                              facecolor="#fca5a5", edgecolor="#dc2626", alpha=0.25, linewidth=1.2)
            ax.add_patch(excl)
        pipe = plt.Circle((1000.0 * xc, 1000.0 * yc), 1000.0 * layout["tube_radius_m"],
                          facecolor="#0f766e", edgecolor="#134e4a", alpha=0.95, linewidth=1.2)
        ax.add_patch(pipe)

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Bed width coordinate (mm)")
    ax.set_ylabel("Filled-bed height coordinate (mm)")
    ax.set_title("Aeration-Tube Cross Section, Worm Exclusion Envelope, and Wall-Return Bands")
    ax.grid(True, alpha=0.15)
    fig.text(
        0.02,
        0.05,
        "Green rectangles show the filled bed cross section. Dark circles are the aeration pipes. "
        "Red circles show the union of the pipe radius and the plume-model safe-distance cutoff. "
        "Orange bands show the additional wall-return unsafe regions from the 1D wall-conduction post-processing. "
        "The reported unsafe volume uses the overlap-aware union of the exclusion cylinders and wall-return bands "
        "clipped to the filled-bed volume.",
        fontsize=9,
        wrap=True,
    )
    fig.savefig(output_dir / "habitat_exclusion_cross_section.png", dpi=220)
    plt.close(fig)


def plot_radial_air_cooling_profile(payload: dict, config: dict, output_dir: Path) -> None:
    radial = compute_radial_profile(payload, config)
    if not radial:
        return

    fig, ax = plt.subplots(figsize=(10, 6.2))
    fig.subplots_adjust(left=0.11, right=0.97, top=0.86, bottom=0.20)
    ax.plot(1000.0 * radial["distances_m"], radial["temperatures_C"], linewidth=2.2, color="#155e75")
    ax.axhline(radial["threshold_C"], color="crimson", linestyle="--", linewidth=1.4)
    if np.isfinite(radial["safe_distance_m"]) and radial["safe_distance_m"] <= radial["distances_m"][-1]:
        ax.axvline(1000.0 * radial["safe_distance_m"], color="crimson", linestyle=":", linewidth=1.4)
        ax.scatter([1000.0 * radial["safe_distance_m"]], [radial["threshold_C"]], color="crimson", zorder=3)
    ax.set_xlabel("Distance from aeration tube outer surface (mm)")
    ax.set_ylabel("Released-air temperature (C)")
    ax.set_title(
        "Radial Air-Cooling Profile Around the Active Operating Point\n"
        f"Profile taken from the active {radial['point_label']} point"
    )
    ax.grid(True, alpha=0.25)

    if radial["threshold_state"] == "not_needed":
        threshold_note = f"Outlet air is already at or below {radial['threshold_C']:.1f} C."
    elif radial["threshold_state"] == "reached":
        threshold_note = (
            f"Safe distance to {radial['threshold_C']:.1f} C: "
            f"{1000.0 * radial['safe_distance_m']:.1f} mm from the tube outer surface."
        )
    elif radial["threshold_state"] == "beyond_plot_limit":
        threshold_note = (
            f"{radial['threshold_C']:.1f} C is reached at about "
            f"{1000.0 * radial['safe_distance_m']:.1f} mm, beyond the plotted range."
        )
    else:
        threshold_note = (
            f"The profile asymptote stays above {radial['threshold_C']:.1f} C, so no finite safe distance "
            "exists in this simple model."
        )

    fig.text(
        0.02,
        0.05,
        "Model: the warm air is treated as a finite top-hat plume that starts with radius b0 and spreads as "
        "b(x) = b0 + alpha x, optionally capped. Temperature and vapor content are re-solved together with "
        "m'_rel c_p dT/dx = -(q'_sens + q'_lat), using a local Lewis-relation evaporation term so the latent sink "
        "falls as the plume cools and humidifies. The contacted wet-bed volume grows gradually with distance "
        "instead of assuming the full local wetted area is active immediately. "
        + threshold_note,
        fontsize=9,
        wrap=True,
    )
    fig.savefig(output_dir / "radial_air_cooling_profile.png", dpi=220)
    plt.close(fig)


def duplicate_series_note(series: list[tuple[str, np.ndarray, np.ndarray]]) -> str | None:
    groups: dict[tuple[tuple[float, ...], tuple[float, ...]], list[str]] = {}
    for label, x_vals, y_vals in series:
        key = (tuple(np.round(x_vals, 6)), tuple(np.round(y_vals, 6)))
        groups.setdefault(key, []).append(label)
    overlapping = [labels for labels in groups.values() if len(labels) > 1]
    if not overlapping:
        return None
    return "; ".join(" = ".join(group) for group in overlapping)


def plot_wire_trade_study(
    payload: dict,
    config: dict,
    output_dir: Path,
    plot_points: int,
    sweep_key: str = "wire_diameter_sweep",
    figure_title: str = "Wire Diameter Trade Study",
    output_name: str = "wire_diameter_trade_study.png",
) -> None:
    raw_entries = payload.get(sweep_key, [])
    active_label = str(payload.get("active_configuration", {}).get("label", ""))
    wire = [entry for entry in raw_entries if normalize_label(str(entry.get("label", ""))) != "covered_top"]
    order_map = {"uncovered_top": 0, "covered_top_plus_vent": 1}
    wire = sorted(
        wire,
        key=lambda entry: (
            1 if same_configuration_label(active_label, str(entry.get("label", ""))) else 0,
            order_map.get(normalize_label(str(entry.get("label", ""))), 2),
        ),
    )
    if not wire:
        return

    first_entry = wire[0]
    x_field = str(first_entry.get("xField", "diameter_mm"))
    x_label = str(first_entry.get("xLabel", "Wire diameter"))
    x_unit = str(first_entry.get("xUnit", "mm"))
    x_axis_label = f"{x_label} ({x_unit})" if x_unit else x_label

    fig, axes = plt.subplots(6, 1, figsize=(14, 19))
    fig.subplots_adjust(left=0.14, right=0.98, top=0.93, bottom=0.10, hspace=0.24)
    fig.suptitle(
        f"{figure_title}\n"
        "Standalone Python rendering from JSON data\n"
        "Solid = hard-safe selected point, dashed = best-available reference point",
        fontsize=15,
        fontweight="bold",
    )

    linestyles = ["-", "--", ":"]
    markers = ["o", "s", "^"]
    series_defs = [
        (
            "selectedPower_W",
            "bestAvailablePower_W",
            "Selected/reference total power (W)",
            "Hard-safe selected power (solid) and best-available reference power (dashed)",
        ),
        (
            "selectedFlow_Lpm",
            "bestAvailableFlow_Lpm",
            "Selected/reference flow rate (L/min)",
            "Hard-safe selected flow (solid) and best-available reference flow (dashed)",
        ),
        (
            "selectedCurrent_A",
            "bestAvailableCurrent_A",
            "Selected/reference current (A)",
            "Hard-safe selected current (solid) and best-available reference current (dashed)",
        ),
        (
            "selectedWireMax_C",
            "bestAvailableWireMax_C",
            "Selected/reference max wire temperature (C)",
            "Hard-safe selected wire temperature (solid) and best-available reference wire temperature (dashed)",
        ),
        (
            "selectedColderNode_C",
            "bestAvailableColderNode_C",
            "Selected/reference colder-node temperature (C)",
            "Hard-safe selected colder-node bed temperature (solid) and best-available reference colder-node temperature (dashed)",
        ),
    ]

    missing_labels: list[str] = []
    duplicate_notes: list[str] = []

    for ax, (selected_field, reference_field, ylabel, title) in zip(axes[:5], series_defs):
        drawn_series: list[tuple[str, np.ndarray, np.ndarray]] = []
        for idx, entry in enumerate(wire):
            diam = as_float_array(entry.get("xValues", entry.get(x_field, [])))
            physical_mask = trade_entry_valid_mask(entry, config)
            if physical_mask.size == len(diam):
                diam = np.where(physical_mask, diam, np.nan)
            y_selected = as_float_array(entry.get(selected_field, []))
            if y_selected.size == 0:
                y_selected = np.full_like(diam, np.nan)
            y_reference = as_float_array(entry.get(reference_field, entry.get(selected_field, [])))
            if y_reference.size == 0:
                y_reference = np.full_like(diam, np.nan)
            if physical_mask.size == len(y_selected):
                y_selected = np.where(physical_mask, y_selected, np.nan)
            if physical_mask.size == len(y_reference):
                y_reference = np.where(physical_mask, y_reference, np.nan)
            selection_modes = entry.get("selectionMode", [])
            if isinstance(selection_modes, list) and len(selection_modes) == len(diam):
                selected_mask = np.array([mode in {"feasible", "fallback"} for mode in selection_modes], dtype=bool)
                y_selected = np.where(selected_mask, y_selected, np.nan)

            valid_idx = np.where(np.isfinite(y_reference))[0]
            if valid_idx.size == 0:
                if entry["label"] not in missing_labels:
                    missing_labels.append(entry["label"])
                continue

            color = np.array(entry["color"], dtype=float)
            is_active = same_configuration_label(active_label, str(entry.get("label", "")))
            linestyle = linestyles[idx % len(linestyles)]
            marker = markers[idx % len(markers)]
            display_label = display_configuration_label(str(entry["label"]))
            ax.plot(
                diam,
                y_reference,
                color=color,
                linewidth=1.4 if not is_active else 1.8,
                linestyle="--",
                alpha=0.75 if not is_active else 0.90,
                zorder=2 + idx,
            )
            ax.plot(
                diam,
                y_selected,
                color=color,
                linewidth=1.8 if not is_active else 2.6,
                linestyle=linestyle,
                label=display_label,
                zorder=5 + idx if is_active else 3 + idx,
            )

            if isinstance(selection_modes, list) and len(selection_modes) == len(diam):
                feasible_idx = np.array([i for i, mode in enumerate(selection_modes) if mode == "feasible"], dtype=int)
                fallback_idx = np.array([i for i, mode in enumerate(selection_modes) if mode == "fallback"], dtype=int)
                if feasible_idx.size > 0:
                    ax.plot(
                        diam[feasible_idx],
                        y_selected[feasible_idx],
                        linestyle="None",
                        marker=marker,
                        markersize=4.8,
                        markerfacecolor="white",
                        markeredgecolor=color,
                        color=color,
                    )
                if fallback_idx.size > 0:
                    ax.plot(
                        diam[fallback_idx],
                        y_selected[fallback_idx],
                        linestyle="None",
                        marker="s",
                        markersize=4.8,
                        markerfacecolor="white",
                        markeredgecolor=color,
                        color=color,
                    )

            best_index = choose_trade_display_index(entry, config, selected=True)
            if best_index is not None:
                if 0 <= best_index < len(diam) and np.isfinite(y_selected[best_index]):
                    ax.plot(
                        diam[best_index],
                        float(y_selected[best_index]),
                        marker="*",
                        markersize=12,
                        color="yellow",
                        markeredgecolor="black",
                    )
            best_available_index = choose_trade_display_index(entry, config, selected=False)
            if best_available_index is not None:
                if 0 <= best_available_index < len(diam) and np.isfinite(y_reference[best_available_index]):
                    ax.plot(
                        diam[best_available_index],
                        float(y_reference[best_available_index]),
                        marker="D",
                        markersize=5.5,
                        color=color,
                        markerfacecolor="white",
                        markeredgecolor="black",
                    )

            drawn_series.append((display_label, diam[np.isfinite(y_reference)], y_reference[np.isfinite(y_reference)]))

        note = duplicate_series_note(drawn_series)
        if note:
            duplicate_notes.append(note)

        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, alpha=0.30)

    ax = axes[5]
    for idx, entry in enumerate(wire):
        diam = as_float_array(entry.get("xValues", entry.get(x_field, [])))
        physical_mask = trade_entry_valid_mask(entry, config)
        if physical_mask.size == len(diam):
            diam = np.where(physical_mask, diam, np.nan)
        color = np.array(entry["color"], dtype=float)
        is_active = same_configuration_label(active_label, str(entry.get("label", "")))
        y_base = float(len(wire) - idx)
        selection_modes = entry.get("selectionMode", [])
        selected_mask = np.isfinite(as_float_array(entry.get("selectedPower_W", [])))
        reference_mask = np.isfinite(as_float_array(entry.get("bestAvailablePower_W", entry.get("selectedPower_W", []))))
        if physical_mask.size == len(selected_mask):
            selected_mask &= physical_mask
        if physical_mask.size == len(reference_mask):
            reference_mask &= physical_mask
        if isinstance(selection_modes, list) and len(selection_modes) == len(diam):
            feasible_mask = np.array([mode == "feasible" for mode in selection_modes], dtype=bool)
            fallback_mask = np.array([mode == "fallback" for mode in selection_modes], dtype=bool)
            feasible_mask &= reference_mask
            fallback_mask &= reference_mask
            missing_selected_mask = reference_mask & ~selected_mask
            if np.any(missing_selected_mask):
                ax.plot(
                    diam[missing_selected_mask],
                    np.full(np.count_nonzero(missing_selected_mask), y_base),
                    linestyle="None",
                    marker="x",
                    markersize=4.0,
                    color=color,
                    alpha=0.75 if not is_active else 0.95,
                )
            if np.any(fallback_mask):
                ax.plot(
                    diam[fallback_mask],
                    np.full(np.count_nonzero(fallback_mask), y_base),
                    linestyle="None",
                    marker="s",
                    markersize=4.4,
                    color=color,
                    markerfacecolor="white",
                    alpha=0.85 if not is_active else 1.0,
                )
            if np.any(feasible_mask):
                ax.plot(
                    diam[feasible_mask],
                    np.full(np.count_nonzero(feasible_mask), y_base),
                    linestyle="None",
                    marker="o",
                    markersize=4.6,
                    color=color,
                    markerfacecolor="white",
                    alpha=0.9 if not is_active else 1.0,
                )
        else:
            feasible = np.array(entry.get("isFeasible", []), dtype=bool)
            ax.plot(
                diam[~feasible],
                np.full(np.count_nonzero(~feasible), y_base),
                linestyle="None",
                marker="x",
                markersize=4,
                color=color,
                alpha=0.65 if not is_active else 0.90,
            )
            ax.plot(
                diam[feasible],
                np.full(np.count_nonzero(feasible), y_base),
                linestyle="None",
                marker="o",
                markersize=4.5,
                color=color,
                markerfacecolor="white",
            )
        best_index = choose_trade_display_index(entry, config, selected=True)
        if best_index is not None:
            if 0 <= best_index < len(diam):
                ax.plot(
                    diam[best_index],
                    y_base,
                    marker="*",
                    markersize=12,
                    color="yellow",
                    markeredgecolor="black",
                )
        best_available_index = choose_trade_display_index(entry, config, selected=False)
        if best_available_index is not None:
            if 0 <= best_available_index < len(diam):
                ax.plot(
                    diam[best_available_index],
                    y_base,
                    marker="D",
                    markersize=5.5,
                    color=color,
                    markerfacecolor="white",
                    markeredgecolor="black",
                )

    ax.set_yticks(
        [float(len(wire) - idx) for idx in range(len(wire))],
        [display_configuration_label(str(entry["label"])) for entry in wire],
    )
    ax.set_xlabel(x_axis_label)
    ax.set_ylabel("Selection status")
    ax.set_title(f"All sampled {x_label.lower()} values (o feasible, s hard-safe fallback, x best-available only, * best selected, D best available)")
    ax.grid(True, alpha=0.25)

    axes[0].legend(loc="best")

    sample_count = len(first_entry.get("xValues", first_entry.get(x_field, [])))
    note_lines = [f"Each series is sampled at {sample_count} {x_label.lower()} values in the exported dataset."] if wire else []
    if missing_labels:
        note_lines.append(
            "No best-available trade-study points were exported for: "
            + ", ".join(sorted(display_configuration_label(label) for label in missing_labels))
        )
    note_lines.append("Solid curves show hard-safe selected points only. Dashed curves show the best-available reference point at each sampled geometry value.")
    note_lines.append("This prevents last-resort points from appearing as if they were valid selected diameters.")
    if active_label:
        note_lines.append(
            "The active configuration is drawn last and slightly thicker so it remains visible when curves overlap."
        )
    if sweep_key == "coil_diameter_sweep":
        heater_inputs = effective_heater_inputs(config.get("summary_inputs", {}), normalized_model_overrides(config))
        wire_inputs = config.get("summary_inputs", {}).get("wire", {})
        max_fit_mm = float(heater_inputs.get("ID_mm", np.nan)) - float(wire_inputs.get("diameter_mm", np.nan))
        if np.isfinite(max_fit_mm):
            for ax in axes:
                ax.axvline(max_fit_mm, color="black", linestyle=":", linewidth=1.1, alpha=0.7)
            note_lines.append(
                f"Physical fit limit for coil mean diameter = heater ID - wire diameter = {max_fit_mm:.3f} mm."
            )
    if duplicate_notes:
        note_lines.append("Overlapping series in this data: " + "; ".join(sorted(set(duplicate_notes))))
    if note_lines:
        fig.text(0.02, 0.05, " | ".join(note_lines), fontsize=9, wrap=True)

    fig.savefig(output_dir / output_name, dpi=220)
    plt.close(fig)


def format_sweep(values: list[float]) -> str:
    arr = np.array(values, dtype=float)
    if arr.size == 0:
        return "n/a"
    if arr.size == 1:
        return f"{arr[0]:.1f}"
    diffs = np.diff(arr)
    step = float(np.median(diffs))
    if np.allclose(diffs, step, rtol=1e-6, atol=1e-9):
        decimals = 0 if np.allclose(arr, np.round(arr)) and np.isclose(step, round(step)) else 1
        return f"{arr[0]:.{decimals}f}:{step:.{decimals}f}:{arr[-1]:.{decimals}f}"
    return f"{arr[0]:.1f}..{arr[-1]:.1f}"


def recommended_hard_safe(rec: dict, limits: dict, opt_cfg: dict) -> bool:
    if not rec:
        return False
    bottom = float(rec.get("TbottomFullPower_C", np.nan))
    top = float(rec.get("TtopFullPower_C", np.nan))
    spread = abs(bottom - top) if np.isfinite(bottom) and np.isfinite(top) else np.inf
    return (
        float(rec["totalCurrent_A"]) <= float(limits["maxTotalCurrent_A"])
        and float(rec["wireMax_C"]) <= float(limits["maxWireTemp_C"])
        and float(rec["airOutlet_C"]) <= float(limits["maxAirOutletTemp_C"])
        and float(rec["deltaP_Pa"]) <= float(limits["maxPressureDrop_Pa"])
        and float(rec.get("waterLoss_kg_day", 0.0)) <= float(limits.get("maxWaterLoss_kg_day", np.inf))
        and float(rec.get("holeVelocity_m_s", 0.0)) <= float(limits.get("maxHoleVelocity_m_s", np.inf))
        and bottom >= float(opt_cfg["min_bottom_temp_C"])
        and top >= float(opt_cfg["min_top_temp_C"])
        and spread <= float(opt_cfg["spread_limit_C"])
        and float(rec["QtoBed_W"]) > 0.0
    )


def recommended_status_short(rec: dict, limits: dict, opt_cfg: dict) -> str:
    if not rec:
        return "none"
    if rec.get("selectionStage") == "feasible" or bool(rec.get("isFeasibleByCriteria", False)):
        return "feasible"
    if rec.get("selectionStage") == "fallback" or bool(rec.get("isConstraintSafeByCriteria", False)):
        return "hard-safe"
    if recommended_hard_safe(rec, limits, opt_cfg):
        return "hard-safe"
    return "best-available only"


def recommended_status_long(rec: dict, limits: dict, opt_cfg: dict) -> str:
    if not rec:
        return "no hard-constraint-safe candidate exists in the sweep"
    if rec.get("selectionStage") == "feasible" or bool(rec.get("isFeasibleByCriteria", False)):
        return "meets all hard constraints and the heating target"
    if rec.get("selectionStage") == "fallback" or bool(rec.get("isConstraintSafeByCriteria", False)):
        return "meets all hard constraints but does not fully meet the heating target"
    if recommended_hard_safe(rec, limits, opt_cfg):
        return "meets all hard constraints but does not fully meet the heating target"
    return "best-available only; it violates one or more hard constraints"


def optimization_config(config: dict) -> dict:
    model_opt = normalized_model_overrides(config).get("optimization", {})
    opt = dict(config.get("optimization", {}))
    opt.setdefault("enabled", True)
    opt.setdefault("top_candidate_count", 8)
    opt.setdefault("n_top", 10)
    opt.setdefault("priority_order", list(DEFAULT_PRIORITY_ORDER))
    opt.setdefault("real_world_guidance", [])
    if "spreadLimit_C" in model_opt:
        opt["spread_limit_C"] = float(model_opt["spreadLimit_C"])
    else:
        opt.setdefault("spread_limit_C", 10.0)
    if "minBottomTemp_C" in model_opt:
        opt["min_bottom_temp_C"] = float(model_opt["minBottomTemp_C"])
    else:
        opt.setdefault("min_bottom_temp_C", 15.0)
    if "minTopTemp_C" in model_opt:
        opt["min_top_temp_C"] = float(model_opt["minTopTemp_C"])
    else:
        opt.setdefault("min_top_temp_C", 15.0)
    return opt


def active_design_points_array(payload: dict) -> list[dict]:
    return payload["design_points"]


def hard_constraint_safe(metrics: dict[str, float | str]) -> bool:
    if "isConstraintSafeByCriteria" in metrics:
        return bool(metrics.get("isConstraintSafeByCriteria", False))
    return all(float(metrics[name]) <= 1e-12 for name in HARD_CONSTRAINT_KEYS) and float(metrics["QtoBed_W"]) > 0.0


def select_recommended_point(payload: dict, config: dict, points: list[dict] | None = None) -> dict:
    ranked = rank_points_by_criteria(payload, config, points)
    hard_safe = [pt for pt in ranked if hard_constraint_safe(pt)]
    return hard_safe[0] if hard_safe else {}


def select_best_available_point(payload: dict, config: dict, points: list[dict] | None = None) -> dict:
    ranked = rank_points_by_criteria(payload, config, points)
    return ranked[0] if ranked else {}


def point_metrics(pt: dict, payload: dict, config: dict) -> dict[str, float | str]:
    limits = effective_limits(payload, config)
    opt_cfg = optimization_config(config)

    bottom = float(pt.get("TbottomFullPower_C", np.nan))
    top = float(pt.get("TtopFullPower_C", np.nan))
    if not np.isfinite(bottom):
        bottom = -1.0e6
    if not np.isfinite(top):
        top = -1.0e6

    spread = abs(bottom - top)
    q_to_bed = float(pt.get("QtoBed_W", 0.0))
    current_excess = max(0.0, float(pt["totalCurrent_A"]) - float(limits["maxTotalCurrent_A"]))
    wire_excess = max(0.0, float(pt["wireMax_C"]) - float(limits["maxWireTemp_C"]))
    air_excess = max(0.0, float(pt["airOutlet_C"]) - float(limits["maxAirOutletTemp_C"]))
    dp_excess = max(0.0, float(pt["deltaP_Pa"]) - float(limits["maxPressureDrop_Pa"]))
    water_loss_excess = max(0.0, float(pt.get("waterLoss_kg_day", 0.0)) - float(limits.get("maxWaterLoss_kg_day", np.inf)))
    hole_velocity = float(pt.get("holeVelocity_m_s", np.nan))
    hole_velocity_excess = max(0.0, hole_velocity - float(limits.get("maxHoleVelocity_m_s", np.inf)))
    heat_shortfall = max(0.0, float(pt["dutyNeeded"]) - 1.0)
    spread_excess = max(0.0, spread - float(opt_cfg["spread_limit_C"]))
    bottom_shortfall = max(0.0, float(opt_cfg["min_bottom_temp_C"]) - bottom)
    top_shortfall = max(0.0, float(opt_cfg["min_top_temp_C"]) - top)
    explicit_safe = (
        current_excess == 0.0
        and wire_excess == 0.0
        and air_excess == 0.0
        and dp_excess == 0.0
        and water_loss_excess == 0.0
        and hole_velocity_excess == 0.0
    )
    thermal_safe = (
        bottom_shortfall == 0.0
        and top_shortfall == 0.0
        and spread_excess == 0.0
        and q_to_bed > 0.0
    )
    constraint_safe = explicit_safe and thermal_safe

    if bool(pt.get("isFeasible", False)) and constraint_safe and heat_shortfall == 0.0:
        stage = "feasible"
        state = "feasible"
    elif constraint_safe:
        stage = "fallback"
        if heat_shortfall > 0.0:
            state = "heat shortfall"
        else:
            state = "constraint-safe"
    elif current_excess > 0.0:
        stage = "last_resort"
        state = "current limit"
    elif wire_excess > 0.0:
        stage = "last_resort"
        state = "wire limit"
    elif air_excess > 0.0:
        stage = "last_resort"
        state = "air limit"
    elif dp_excess > 0.0:
        stage = "last_resort"
        state = "pressure limit"
    elif water_loss_excess > 0.0:
        stage = "last_resort"
        state = "water loss"
    elif hole_velocity_excess > 0.0:
        stage = "last_resort"
        state = "hole velocity"
    elif bottom_shortfall > 0.0:
        stage = "last_resort"
        state = "bottom temperature"
    elif top_shortfall > 0.0:
        stage = "last_resort"
        state = "top temperature"
    elif spread_excess > 0.0:
        stage = "last_resort"
        state = "temperature spread"
    else:
        stage = "last_resort"
        state = "nonpositive heat" if q_to_bed <= 0.0 else "heat shortfall"

    return {
        "current_excess_A": current_excess,
        "wire_temp_excess_C": wire_excess,
        "air_outlet_excess_C": air_excess,
        "pressure_drop_excess_Pa": dp_excess,
        "water_loss_excess_kg_day": water_loss_excess,
        "hole_velocity_excess_m_s": hole_velocity_excess,
        "heat_shortfall": heat_shortfall,
        "temp_spread_excess_C": spread_excess,
        "bottom_temp_shortfall_C": bottom_shortfall,
        "top_temp_shortfall_C": top_shortfall,
        "minus_min_bed_temp_C": -min(bottom, top),
        "minus_mean_bed_temp_C": -0.5 * (bottom + top),
        "power_W": float(pt["totalPower_W"]),
        "temp_spread_C": spread,
        "flowPerHole_Lpm": float(pt.get("flowPerHole_Lpm", np.nan)),
        "holeVelocity_m_s": float(pt.get("holeVelocity_m_s", np.nan)),
        "waterLoss_kg_day": float(pt.get("waterLoss_kg_day", np.nan)),
        "latentEvap_W": float(pt.get("latentEvap_W", np.nan)),
        "QtoBed_W": q_to_bed,
        "bottom_temp_C": bottom,
        "top_temp_C": top,
        "colder_node_temp_C": min(bottom, top),
        "mean_bed_temp_C": 0.5 * (bottom + top),
        "selectionStage": stage,
        "constraintState": state,
        "isConstraintSafeByCriteria": constraint_safe,
        "isFeasibleByCriteria": stage == "feasible",
    }


def criteria_sort_key(metrics: dict[str, float | str], config: dict) -> tuple[float, ...]:
    order = optimization_config(config)["priority_order"]
    return tuple(float(metrics[name]) for name in order)


def rank_points_by_criteria(payload: dict, config: dict, points: list[dict] | None = None) -> list[dict]:
    source = active_design_points_array(payload) if points is None else points
    ranked: list[dict] = []
    for pt in source:
        metrics = point_metrics(pt, payload, config)
        enriched = dict(pt)
        enriched["criteriaKey"] = criteria_sort_key(metrics, config)
        enriched["selectionStage"] = str(metrics["selectionStage"])
        enriched["constraintState"] = str(metrics["constraintState"])
        enriched["isConstraintSafeByCriteria"] = bool(metrics["isConstraintSafeByCriteria"])
        enriched["isFeasibleByCriteria"] = bool(metrics["isFeasibleByCriteria"])
        enriched["tempSpread_C"] = float(metrics["temp_spread_C"])
        enriched["TbottomFullPower_C"] = float(metrics["bottom_temp_C"])
        enriched["TtopFullPower_C"] = float(metrics["top_temp_C"])
        enriched["colderNodeTemp_C"] = float(metrics["colder_node_temp_C"])
        enriched["meanBedTemp_C"] = float(metrics["mean_bed_temp_C"])
        enriched["heatShortfall"] = float(metrics["heat_shortfall"])
        enriched["tempSpreadExcess_C"] = float(metrics["temp_spread_excess_C"])
        enriched["bottomTempShortfall_C"] = float(metrics["bottom_temp_shortfall_C"])
        enriched["topTempShortfall_C"] = float(metrics["top_temp_shortfall_C"])
        enriched["flowPerHole_Lpm"] = float(metrics["flowPerHole_Lpm"])
        enriched["holeVelocity_m_s"] = float(metrics["holeVelocity_m_s"])
        enriched["waterLoss_kg_day"] = float(metrics["waterLoss_kg_day"])
        enriched["latentEvap_W"] = float(metrics["latentEvap_W"])
        ranked.append(enriched)
    ranked.sort(key=lambda pt: pt["criteriaKey"])
    return ranked


def optimization_ranked_points(payload: dict, config: dict) -> list[dict]:
    opt_cfg = optimization_config(config)
    if not opt_cfg.get("enabled", True):
        return []
    ranked = rank_points_by_criteria(payload, config)
    return ranked[: int(opt_cfg.get("n_top", 10))]


def matlab_top_candidates(payload: dict, config: dict, n_top: int | None = None) -> list[dict]:
    opt_cfg = optimization_config(config)
    limit = int(opt_cfg.get("top_candidate_count", 8) if n_top is None else n_top)
    ranked = rank_points_by_criteria(payload, config)
    return ranked[: min(limit, len(ranked))]


def cooling_flow_array(cooling_cfg: dict) -> np.ndarray:
    direct = cooling_cfg.get("totalFlow_Lpm")
    if isinstance(direct, list) and direct:
        return np.array([float(x) for x in direct], dtype=float)
    count = max(2, int(cooling_cfg.get("totalFlow_Lpm_count", 100)))
    return np.linspace(
        float(cooling_cfg.get("totalFlow_Lpm_min", 40.0)),
        float(cooling_cfg.get("totalFlow_Lpm_max", 600.0)),
        count,
    )


def cfm_to_m3_s(flow_cfm: float) -> float:
    return float(flow_cfm) * 0.00047194745


def cfm_to_lpm(flow_cfm: float) -> float:
    return float(flow_cfm) * 28.316846592


def lpm_to_cfm(flow_lpm: float) -> float:
    return float(flow_lpm) / 28.316846592


def in_h2o_to_pa(delta_p_in_h2o: float) -> float:
    return float(delta_p_in_h2o) * 249.08891


def assist_blower_state(cooling_cfg: dict, total_flow_m3_s: float) -> dict:
    blower = dict(cooling_cfg.get("assist_blower", {}))
    enabled = bool(blower.get("enabled", False))
    flow_cfm = total_flow_m3_s / 0.00047194745
    rated_flow_cfm = float(blower.get("ratedFlow_CFM", np.nan))
    rated_pressure_pa = in_h2o_to_pa(float(blower.get("ratedPressure_inH2O", np.nan)))
    shutoff_pressure_pa = in_h2o_to_pa(float(blower.get("shutoffPressure_inH2O", np.nan)))
    free_air_cfm_raw = blower.get("freeAirFlow_CFM", None)
    free_air_cfm = float(free_air_cfm_raw) if free_air_cfm_raw not in (None, "") else float("nan")
    if not np.isfinite(free_air_cfm) and np.isfinite(rated_flow_cfm) and np.isfinite(rated_pressure_pa) and np.isfinite(shutoff_pressure_pa):
        if shutoff_pressure_pa > rated_pressure_pa and rated_flow_cfm > 0.0:
            free_air_cfm = rated_flow_cfm * shutoff_pressure_pa / max(shutoff_pressure_pa - rated_pressure_pa, 1.0e-12)
        else:
            free_air_cfm = rated_flow_cfm
    available_pressure_pa = float("nan")
    if enabled and np.isfinite(free_air_cfm) and free_air_cfm > 0.0 and np.isfinite(shutoff_pressure_pa):
        available_pressure_pa = max(shutoff_pressure_pa * (1.0 - flow_cfm / free_air_cfm), 0.0)
    motor_power_W = float(blower.get("motorPower_W", np.nan))
    electrical_power_W = motor_power_W if enabled and flow_cfm > 0.0 and np.isfinite(motor_power_W) else 0.0
    return {
        "enabled": enabled,
        "model": str(blower.get("model", "")),
        "referenceNote": str(blower.get("reference", "")),
        "ratedFlow_CFM": float(rated_flow_cfm) if np.isfinite(rated_flow_cfm) else float("nan"),
        "ratedPressure_Pa": float(rated_pressure_pa) if np.isfinite(rated_pressure_pa) else float("nan"),
        "shutoffPressure_Pa": float(shutoff_pressure_pa) if np.isfinite(shutoff_pressure_pa) else float("nan"),
        "impliedFreeAirFlow_CFM": float(free_air_cfm) if np.isfinite(free_air_cfm) else float("nan"),
        "requestedFlow_CFM": float(flow_cfm),
        "availablePressure_Pa": float(available_pressure_pa) if np.isfinite(available_pressure_pa) else float("nan"),
        "electricalPower_W": float(electrical_power_W),
    }


def summer_evaporation_loss_segment(
    dm_release_kg_s: float,
    tair_C: float,
    tbed_C: float,
    q_avail_W: float,
    hole_velocity_m_s: float,
    seg_area_m2: float,
    pressure_Pa: float,
    inlet_relative_humidity: float,
    surface_temp_C: float,
    lewis_factor: float,
    evaporation_multiplier: float,
    hole_diameter_m: float,
    inlet_air_ref_C: float,
) -> tuple[float, float]:
    lower_T_C = min(tair_C, tbed_C)
    upper_T_C = max(tair_C, tbed_C)
    surfaceT_C = min(max(surface_temp_C, lower_T_C), upper_T_C)
    surfaceT_K = surfaceT_C + 273.15
    props = air_props(tair_C, pressure_Pa)
    rho_v_surface = water_vapor_density_at_saturation_kg_m3(surfaceT_K)
    rho_v_inlet = inlet_relative_humidity * water_vapor_density_at_saturation_kg_m3(inlet_air_ref_C + 273.15)
    delta_rho_v = max(rho_v_surface - rho_v_inlet, 0.0)

    re_hole = props["rho_kg_m3"] * max(hole_velocity_m_s, 0.0) * max(hole_diameter_m, 1e-12) / max(props["mu_Pa_s"], 1e-12)
    nu_hole = churchill_bernstein_nu(re_hole, props["Pr"])
    h_evap = evaporation_multiplier * nu_hole * props["k_W_mK"] / max(hole_diameter_m, 1e-12)
    h_mass = h_evap / max(props["rho_kg_m3"] * props["cp_J_kgK"] * lewis_factor ** (2.0 / 3.0), 1e-12)

    vol_release_m3_s = dm_release_kg_s / max(props["rho_kg_m3"], 1e-12)
    hfg = latent_heat_vaporization_water_J_kg(surfaceT_K)
    m_diff = h_mass * seg_area_m2 * delta_rho_v
    m_capacity = vol_release_m3_s * delta_rho_v
    m_heat = q_avail_W / max(hfg, 1e-12)
    m_evap = max(0.0, min(m_diff, m_capacity, m_heat))
    return m_evap, m_evap * hfg


def spot_cooler_supply_state(
    cooling_cfg: dict,
    total_flow_m3_s: float,
    ambient_air_C: float,
    pressure_Pa: float,
) -> dict:
    spot = dict(cooling_cfg.get("spot_cooler", {}))
    ambient_rh = float(cooling_cfg.get("inletRelativeHumidity", 0.5))
    air_in = air_props(ambient_air_C, pressure_Pa)
    mdot_total_kg_s = air_in["rho_kg_m3"] * total_flow_m3_s
    cp_air = float(air_in["cp_J_kgK"])
    derive_from_specs = bool(spot.get("deriveTargetFromRatedSpecs", False))
    target_raw = spot.get("targetSupplyAir_C", None)
    target_supply_C = float(target_raw) if target_raw not in (None, "") else ambient_air_C

    rated_total_btu_h_raw = spot.get("ratedTotalCapacity_BTU_h", None)
    rated_total_capacity_W = (
        float(rated_total_btu_h_raw) * 0.29307107 if rated_total_btu_h_raw not in (None, "") else float("nan")
    )
    rated_airflow_cfm_raw = spot.get("ratedAirflow_CFM", None)
    rated_airflow_cfm = float(rated_airflow_cfm_raw) if rated_airflow_cfm_raw not in (None, "") else float("nan")
    rated_airflow_m3_s = rated_airflow_cfm * 0.00047194745 if np.isfinite(rated_airflow_cfm) else float("nan")
    rated_dehum_raw = spot.get("ratedDehumidification_pints_day", None)
    rated_condensate_kg_s = (
        float(rated_dehum_raw) * 0.473176473 / 86400.0 if rated_dehum_raw not in (None, "") else float("nan")
    )
    hfg_ambient = latent_heat_vaporization_water_J_kg(ambient_air_C + 273.15)
    rated_latent_W = rated_condensate_kg_s * hfg_ambient if np.isfinite(rated_condensate_kg_s) else float("nan")
    rated_sensible_W = (
        max(rated_total_capacity_W - rated_latent_W, 0.0)
        if np.isfinite(rated_total_capacity_W) and np.isfinite(rated_latent_W)
        else float("nan")
    )

    derived_target_supply_C = float("nan")
    if np.isfinite(rated_sensible_W) and np.isfinite(rated_airflow_m3_s) and rated_airflow_m3_s > 0.0:
        rated_mdot_kg_s = air_in["rho_kg_m3"] * rated_airflow_m3_s
        if rated_mdot_kg_s > 0.0 and cp_air > 0.0:
            derived_target_supply_C = ambient_air_C - rated_sensible_W / (rated_mdot_kg_s * cp_air)
    if derive_from_specs and np.isfinite(derived_target_supply_C):
        target_supply_C = float(derived_target_supply_C)

    requested_load_W = max(mdot_total_kg_s * cp_air * max(ambient_air_C - target_supply_C, 0.0), 0.0)

    max_capacity_raw = spot.get("maxCoolingCapacity_W", None)
    if max_capacity_raw not in (None, ""):
        max_capacity_W = float(max_capacity_raw)
    elif np.isfinite(rated_sensible_W):
        max_capacity_W = float(rated_sensible_W)
    else:
        max_capacity_W = np.inf
    actual_load_W = min(requested_load_W, max_capacity_W)
    if mdot_total_kg_s > 0.0 and cp_air > 0.0:
        actual_supply_C = ambient_air_C - actual_load_W / (mdot_total_kg_s * cp_air)
    else:
        actual_supply_C = ambient_air_C

    inlet_w = humidity_ratio_from_rh(ambient_air_C, ambient_rh, pressure_Pa)
    inlet_dewpoint_C = dew_point_from_rh(ambient_air_C, ambient_rh, pressure_Pa)
    mdot_dry_air_kg_s = mdot_total_kg_s / max(1.0 + inlet_w, 1e-12)
    saturated_outlet_w = humidity_ratio_saturated(actual_supply_C, pressure_Pa)
    condensate_psych_kg_s = max(mdot_dry_air_kg_s * max(inlet_w - saturated_outlet_w, 0.0), 0.0)
    if np.isfinite(rated_condensate_kg_s):
        condensate_kg_s = min(condensate_psych_kg_s, rated_condensate_kg_s)
    else:
        condensate_kg_s = condensate_psych_kg_s
    outlet_w = max(inlet_w - condensate_kg_s / max(mdot_dry_air_kg_s, 1e-12), 0.0)
    supply_rh = relative_humidity_from_humidity_ratio(actual_supply_C, outlet_w, pressure_Pa)

    cop_raw = spot.get("COP", None)
    if cop_raw in (None, ""):
        power_W = float("nan")
    else:
        cop = float(cop_raw)
        power_W = actual_load_W / cop if cop > 0.0 else float("nan")

    return {
        "targetSupplyAir_C": target_supply_C,
        "derivedTargetSupply_C": float(derived_target_supply_C),
        "actualSupplyAir_C": float(actual_supply_C),
        "coolingLoad_W": float(actual_load_W),
        "coolingPower_W": float(power_W),
        "requestedCoolingLoad_W": float(requested_load_W),
        "massFlow_kg_s": float(mdot_total_kg_s),
        "cpAir_J_kgK": cp_air,
        "ambientRelativeHumidity": float(ambient_rh),
        "supplyRelativeHumidity": float(supply_rh),
        "condensate_kg_day": float(condensate_kg_s * 86400.0),
        "psychrometricCondensate_kg_day": float(condensate_psych_kg_s * 86400.0),
        "ratedCondensate_kg_day": float(rated_condensate_kg_s * 86400.0) if np.isfinite(rated_condensate_kg_s) else float("nan"),
        "ratedTotalCapacity_W": float(rated_total_capacity_W) if np.isfinite(rated_total_capacity_W) else float("nan"),
        "ratedLatentCapacity_W": float(rated_latent_W) if np.isfinite(rated_latent_W) else float("nan"),
        "ratedSensibleCapacity_W": float(rated_sensible_W) if np.isfinite(rated_sensible_W) else float("nan"),
        "ratedAirflow_CFM": float(rated_airflow_cfm) if np.isfinite(rated_airflow_cfm) else float("nan"),
        "ratedDehumidification_pints_day": float(rated_dehum_raw) if rated_dehum_raw not in (None, "") else float("nan"),
        "deriveTargetFromRatedSpecs": derive_from_specs,
        "capacityBasis": str(spot.get("capacityBasis", "")),
        "referenceNote": str(spot.get("reference", "")),
        "ambientDewPoint_C": float(inlet_dewpoint_C),
    }


def simulate_summer_cooling_point(
    config: dict,
    flow_Lpm: float,
    bin_model: dict,
    required_cooling_W: float,
    ambient_air_override_C: float | None = None,
) -> dict:
    summary_inputs = config.get("summary_inputs", {})
    model_overrides = normalized_model_overrides(config)
    cooling_cfg = cooling_mode_config(config)
    env_inputs = effective_environment_inputs(summary_inputs, model_overrides)
    ambient_air_C = float(cooling_cfg["ambientAir_C"] if ambient_air_override_C is None else ambient_air_override_C)
    pressure_Pa = float(cooling_cfg.get("pressure_Pa", env_inputs.get("pressure_Pa", 101325.0)))
    design_bed_C = float(cooling_cfg["designBed_C"])
    bin_inputs = effective_bin_inputs(summary_inputs, model_overrides)
    aer_inputs = effective_aeration_inputs(summary_inputs, model_overrides)
    wetted_area_m2 = derived_wetted_area_m2(summary_inputs, model_overrides)
    evaporation_inputs = model_overrides.get("evaporation", {})
    calibration = model_overrides.get("calibration", {})

    n_tubes = max(1, int(round(float(aer_inputs.get("nParallelTubes", 1)))))
    aer_len_m = max(1e-12, float(aer_inputs.get("length_m", 1.22)))
    aer_id_m = max(1e-12, float(aer_inputs.get("ID_mm", 40.0)) / 1000.0)
    aer_od_m = max(aer_id_m, float(aer_inputs.get("OD_mm", 44.0)) / 1000.0)
    wall_k = float(aer_inputs.get("wallK_W_mK", 16.0))
    segments = max(1, int(round(float(aer_inputs.get("segments", 40)))))
    bed_h = float(aer_inputs.get("bedH_W_m2K", 8.0)) * float(calibration.get("bedHTMultiplier", 1.0))
    release_fraction = float(aer_inputs.get("releaseFraction", 0.95))
    end_minor_k = float(aer_inputs.get("endMinorK", 1.0))

    hole_diameter_m = max(1e-12, float(model_overrides.get("perforation", {}).get("holeDiameter_m", 0.005)))
    holes_per_tube = max(1, int(round(float(model_overrides.get("perforation", {}).get("holesPerTube", 60)))))
    ambient_relative_humidity = float(cooling_cfg.get("inletRelativeHumidity", evaporation_inputs.get("relativeHumidity", 0.80)))
    surface_temp_C = float(cooling_cfg.get("surfaceTemp_C", evaporation_inputs.get("surfaceTemp_C", design_bed_C)))
    lewis_factor = float(evaporation_inputs.get("lewisFactor", 1.0))
    evaporation_multiplier = float(calibration.get("evaporationHMultiplier", 1.0))
    bottom_fraction = float(bin_inputs.get("heatToBottomFraction", 0.60))

    total_flow_m3_s = flow_Lpm / 60000.0
    air_in = air_props(ambient_air_C, pressure_Pa)
    mdot_tube_kg_s = air_in["rho_kg_m3"] * total_flow_m3_s / n_tubes
    spot = spot_cooler_supply_state(cooling_cfg, total_flow_m3_s, ambient_air_C, pressure_Pa)
    assist = assist_blower_state(cooling_cfg, total_flow_m3_s)
    supply_air_C = float(spot["actualSupplyAir_C"])
    supply_relative_humidity = float(spot.get("supplyRelativeHumidity", ambient_relative_humidity))
    supply_props = air_props(supply_air_C, pressure_Pa)
    atube_m2 = np.pi * aer_id_m**2 / 4.0
    pout_m = np.pi * aer_od_m
    dx_m = aer_len_m / segments
    released_total = mdot_tube_kg_s * release_fraction
    hole_area_m2 = np.pi * hole_diameter_m**2 / 4.0
    holes_per_segment = max(float(holes_per_tube) / max(segments, 1), 1e-12)
    vol_release_m3_s = released_total / max(supply_props["rho_kg_m3"], 1e-12)
    flow_per_hole_lpm = vol_release_m3_s * 60000.0 / holes_per_tube
    hole_velocity_m_s = vol_release_m3_s / max(holes_per_tube * hole_area_m2, 1e-12)
    seg_area_m2 = wetted_area_m2 / max(n_tubes * segments, 1)

    def branch_response(bed_ref_C: float) -> dict:
        remaining_mdot = mdot_tube_kg_s
        dm_release = released_total / segments
        tair_C = supply_air_C
        q_to_bed = 0.0
        latent_evap = 0.0
        water_loss_kg_s = 0.0
        dp_tube = 0.0
        local_hole_velocity = hole_velocity_m_s
        air_profile = np.zeros(segments + 1, dtype=float)
        air_profile[0] = tair_C

        for idx in range(segments):
            props = air_props(tair_C, pressure_Pa)
            u_bulk = remaining_mdot / max(props["rho_kg_m3"] * atube_m2, 1e-12)
            re_tube = props["rho_kg_m3"] * u_bulk * aer_id_m / max(props["mu_Pa_s"], 1e-12)
            nu_tube = internal_tube_nu(re_tube, props["Pr"], aer_id_m, aer_len_m)
            h_inside = nu_tube * props["k_W_mK"] / aer_id_m
            u_tube = 1.0 / max(
                1.0 / max(h_inside, 1e-12)
                + (aer_od_m - aer_id_m) / (2.0 * max(wall_k, 1e-12))
                + 1.0 / max(bed_h, 1e-12),
                1e-12,
            )
            dm_kg_s = min(dm_release, remaining_mdot)
            remaining_after = max(remaining_mdot - dm_kg_s, 1e-12)
            mdot_mean = max(remaining_mdot - 0.5 * dm_kg_s, 1e-12)
            ua_seg_W_K = u_tube * pout_m * dx_m
            wall_factor = np.exp(-ua_seg_W_K / max(mdot_mean * props["cp_J_kgK"], 1e-12))
            tair_after_wall_C = bed_ref_C + (tair_C - bed_ref_C) * wall_factor
            q_wall = mdot_mean * props["cp_J_kgK"] * (tair_C - tair_after_wall_C)
            t_release_C = 0.5 * (tair_C + tair_after_wall_C)
            q_jet = dm_kg_s * props["cp_J_kgK"] * (t_release_C - bed_ref_C)
            release_props = air_props(t_release_C, pressure_Pa)
            vol_release_seg_m3_s = dm_kg_s / max(release_props["rho_kg_m3"], 1e-12)
            local_hole_velocity = vol_release_seg_m3_s / max(holes_per_segment * hole_area_m2, 1e-12)
            q_avail = dm_kg_s * release_props["cp_J_kgK"] * abs(t_release_C - bed_ref_C)
            m_evap, q_evap = summer_evaporation_loss_segment(
                dm_kg_s,
                t_release_C,
                bed_ref_C,
                q_avail,
                local_hole_velocity,
                seg_area_m2,
                pressure_Pa,
                supply_relative_humidity,
                surface_temp_C,
                lewis_factor,
                evaporation_multiplier,
                hole_diameter_m,
                supply_air_C,
            )
            q_to_bed += q_wall + q_jet - q_evap
            latent_evap += q_evap
            water_loss_kg_s += m_evap
            tair_C = tair_after_wall_C
            remaining_mdot = remaining_after
            air_profile[idx + 1] = tair_C
            f = churchill_friction_factor(re_tube)
            dp_tube += (f * dx_m / aer_id_m + end_minor_k / segments) * 0.5 * props["rho_kg_m3"] * u_bulk**2

        distribution_dp = branch_distribution_losses(
            total_flow_m3_s,
            total_flow_m3_s / max(n_tubes, 1),
            supply_air_C,
            pressure_Pa,
            n_tubes,
            aer_inputs,
            aer_id_m,
        )
        return {
            "q_to_bed_W": float(q_to_bed),
            "latent_evap_W": float(latent_evap),
            "water_loss_kg_s": float(water_loss_kg_s),
            "delta_p_Pa": float(dp_tube + distribution_dp["total_Pa"]),
            "air_profile_C": air_profile,
            "air_outlet_C": float(tair_C),
            "hole_velocity_m_s": float(local_hole_velocity),
            "distributionDeltaP_Pa": float(distribution_dp["total_Pa"]),
            "distributionHeader_Pa": float(distribution_dp["header_Pa"]),
            "distributionHeaderToSplitter_Pa": float(distribution_dp["header_to_splitter_Pa"]),
            "distributionSplitterBody_Pa": float(distribution_dp["splitter_body_Pa"]),
            "distributionConnectorFriction_Pa": float(distribution_dp["connector_friction_Pa"]),
            "distributionConnectorToBranch_Pa": float(distribution_dp["connector_to_branch_Pa"]),
        }

    lower_bed_ref_C = min(design_bed_C, supply_air_C, ambient_air_C)
    upper_bed_ref_C = max(design_bed_C, supply_air_C, ambient_air_C)

    def bed_reference_residual(bed_ref_trial_C: float) -> tuple[float, dict, float, float]:
        trial_state = branch_response(bed_ref_trial_C)
        tb_trial_C, tt_trial_C = two_node_steady_state_python(
            bin_model, ambient_air_C, trial_state["q_to_bed_W"], bottom_fraction
        )
        mean_trial_C = 0.5 * (tb_trial_C + tt_trial_C)
        return mean_trial_C - bed_ref_trial_C, trial_state, tb_trial_C, tt_trial_C

    res_lo, state_lo, tb_lo, tt_lo = bed_reference_residual(lower_bed_ref_C)
    res_hi, state_hi, tb_hi, tt_hi = bed_reference_residual(upper_bed_ref_C)

    if np.sign(res_lo) == np.sign(res_hi):
        candidates = [
            (abs(res_lo), lower_bed_ref_C, state_lo, tb_lo, tt_lo),
            (abs(res_hi), upper_bed_ref_C, state_hi, tb_hi, tt_hi),
        ]
        _, bed_ref_C, state, tb_C, tt_C = min(candidates, key=lambda item: item[0])
    else:
        state = state_lo
        tb_C = tb_lo
        tt_C = tt_lo
        lo = lower_bed_ref_C
        hi = upper_bed_ref_C
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            res_mid, state_mid, tb_mid, tt_mid = bed_reference_residual(mid)
            bed_ref_C = mid
            state = state_mid
            tb_C = tb_mid
            tt_C = tt_mid
            if abs(res_mid) < 0.02 or abs(hi - lo) < 0.02:
                break
            if np.sign(res_mid) == np.sign(res_lo):
                lo = mid
                res_lo = res_mid
            else:
                hi = mid
                res_hi = res_mid
    q_to_bed = float(state["q_to_bed_W"])
    latent_evap = float(state["latent_evap_W"])
    water_loss_kg_s = float(state["water_loss_kg_s"])
    delta_p = float(state["delta_p_Pa"])
    tair_C = float(state["air_outlet_C"])
    air_profile = state["air_profile_C"]
    hole_velocity_m_s = float(state["hole_velocity_m_s"])

    cooling_capacity_W = max(-q_to_bed, 0.0)
    cooling_shortfall = max(required_cooling_W - cooling_capacity_W, 0.0)
    cooling_duty = required_cooling_W / max(cooling_capacity_W, 1e-12) if cooling_capacity_W > 0 else np.inf

    return {
        "pitch_mm": 0.0,
        "voltage_V": 0.0,
        "totalFlow_Lpm": float(flow_Lpm),
        "totalPower_W": 0.0,
        "totalCurrent_A": 0.0,
        "wireMax_C": ambient_air_C,
        "airOutlet_C": float(tair_C),
        "airProfile_C": air_profile.tolist(),
        "airInlet_C": float(supply_air_C),
        "airInletRelativeHumidity": float(supply_relative_humidity),
        "ambientRelativeHumidity": float(ambient_relative_humidity),
        "deltaP_Pa": float(delta_p),
        "distributionDeltaP_Pa": float(state.get("distributionDeltaP_Pa", np.nan)),
        "distributionHeader_Pa": float(state.get("distributionHeader_Pa", np.nan)),
        "distributionHeaderToSplitter_Pa": float(state.get("distributionHeaderToSplitter_Pa", np.nan)),
        "distributionSplitterBody_Pa": float(state.get("distributionSplitterBody_Pa", np.nan)),
        "distributionConnectorFriction_Pa": float(state.get("distributionConnectorFriction_Pa", np.nan)),
        "distributionConnectorToBranch_Pa": float(state.get("distributionConnectorToBranch_Pa", np.nan)),
        "flowPerHole_Lpm": float(flow_per_hole_lpm),
        "holeVelocity_m_s": float(hole_velocity_m_s),
        "waterLoss_kg_day": float(water_loss_kg_s * 86400.0),
        "latentEvap_W": float(latent_evap),
        "QtoBed_W": float(q_to_bed),
        "coolingCapacity_W": float(cooling_capacity_W),
        "coolingShortfall_W": float(cooling_shortfall),
        "coolingDutyNeeded": float(cooling_duty),
        "spotCoolerTargetSupply_C": float(spot["targetSupplyAir_C"]),
        "spotCoolerActualSupply_C": float(spot["actualSupplyAir_C"]),
        "spotCoolerLoad_W": float(spot["coolingLoad_W"]),
        "spotCoolerPower_W": float(spot["coolingPower_W"]),
        "assistBlowerAvailablePressure_Pa": float(assist.get("availablePressure_Pa", np.nan)),
        "assistBlowerPower_W": float(assist.get("electricalPower_W", np.nan)),
        "assistBlowerRatedFlow_CFM": float(assist.get("ratedFlow_CFM", np.nan)),
        "assistBlowerRatedPressure_Pa": float(assist.get("ratedPressure_Pa", np.nan)),
        "assistBlowerShutoffPressure_Pa": float(assist.get("shutoffPressure_Pa", np.nan)),
        "assistBlowerImpliedFreeAirFlow_CFM": float(assist.get("impliedFreeAirFlow_CFM", np.nan)),
        "assistBlowerModel": str(assist.get("model", "")),
        "assistBlowerReferenceNote": str(assist.get("referenceNote", "")),
        "spotCoolerRequestedLoad_W": float(spot["requestedCoolingLoad_W"]),
        "spotCoolerDerivedTargetSupply_C": float(spot.get("derivedTargetSupply_C", np.nan)),
        "spotCoolerSupplyRelativeHumidity": float(spot.get("supplyRelativeHumidity", np.nan)),
        "spotCoolerCondensate_kg_day": float(spot.get("condensate_kg_day", np.nan)),
        "spotCoolerPsychrometricCondensate_kg_day": float(spot.get("psychrometricCondensate_kg_day", np.nan)),
        "spotCoolerRatedCondensate_kg_day": float(spot.get("ratedCondensate_kg_day", np.nan)),
        "spotCoolerRatedTotalCapacity_W": float(spot.get("ratedTotalCapacity_W", np.nan)),
        "spotCoolerRatedLatentCapacity_W": float(spot.get("ratedLatentCapacity_W", np.nan)),
        "spotCoolerRatedSensibleCapacity_W": float(spot.get("ratedSensibleCapacity_W", np.nan)),
        "spotCoolerRatedAirflow_CFM": float(spot.get("ratedAirflow_CFM", np.nan)),
        "spotCoolerRatedDehumidification_pints_day": float(spot.get("ratedDehumidification_pints_day", np.nan)),
        "spotCoolerAmbientDewPoint_C": float(spot.get("ambientDewPoint_C", np.nan)),
        "spotCoolerCapacityBasis": str(spot.get("capacityBasis", "")),
        "spotCoolerReferenceNote": str(spot.get("referenceNote", "")),
        "bedReferenceSolved_C": float(bed_ref_C),
        "TbottomFullPower_C": float(tb_C),
        "TtopFullPower_C": float(tt_C),
        "dutyNeeded": float(cooling_duty),
    }


def cooling_point_metrics(pt: dict, config: dict) -> dict[str, float | str]:
    cooling_cfg = cooling_mode_config(config)
    limits = cooling_cfg["limits"]
    opt = cooling_cfg["optimization"]
    bottom = float(pt.get("TbottomFullPower_C", np.nan))
    top = float(pt.get("TtopFullPower_C", np.nan))
    if not np.isfinite(bottom):
        bottom = 1.0e6
    if not np.isfinite(top):
        top = 1.0e6
    spread = abs(bottom - top)
    pressure_excess = max(0.0, float(pt.get("deltaP_Pa", 0.0)) - float(limits.get("maxPressureDrop_Pa", np.inf)))
    assist_blower_pressure_excess = max(
        0.0,
        float(pt.get("deltaP_Pa", 0.0)) - float(pt.get("assistBlowerAvailablePressure_Pa", np.inf)),
    )
    water_loss = float(pt.get("waterLoss_kg_day", 0.0))
    water_excess = max(0.0, water_loss - float(limits.get("maxWaterLoss_kg_day", np.inf)))
    hole_velocity = float(pt.get("holeVelocity_m_s", np.nan))
    hole_velocity_excess = max(0.0, hole_velocity - float(limits.get("maxHoleVelocity_m_s", np.inf)))
    bottom_excess = max(0.0, bottom - float(opt["max_bottom_temp_C"]))
    top_excess = max(0.0, top - float(opt["max_top_temp_C"]))
    spread_excess = max(0.0, spread - float(opt["spread_limit_C"]))
    cooling_capacity = max(0.0, -float(pt.get("QtoBed_W", 0.0)))
    spot_cooler_load = float(pt.get("spotCoolerLoad_W", np.nan))
    spot_cooler_power = float(pt.get("spotCoolerPower_W", np.nan))
    assist_blower_power = float(pt.get("assistBlowerPower_W", np.nan))
    combined_cooling_power = 0.0
    for value in (spot_cooler_power, assist_blower_power):
        if np.isfinite(value):
            combined_cooling_power += value
    hard_safe = (
        assist_blower_pressure_excess == 0.0
        and pressure_excess == 0.0
        and water_excess == 0.0
        and hole_velocity_excess == 0.0
        and bottom_excess == 0.0
        and top_excess == 0.0
        and spread_excess == 0.0
        and cooling_capacity > 0.0
    )
    if hard_safe:
        stage = "feasible"
        state = "cooling-safe"
    elif assist_blower_pressure_excess > 0.0:
        stage = "last_resort"
        state = "assist-blower static limit"
    elif pressure_excess > 0.0:
        stage = "last_resort"
        state = "pressure limit"
    elif water_excess > 0.0:
        stage = "last_resort"
        state = "water loss"
    elif hole_velocity_excess > 0.0:
        stage = "last_resort"
        state = "hole velocity"
    elif bottom_excess > 0.0:
        stage = "last_resort"
        state = "bottom temperature"
    elif top_excess > 0.0:
        stage = "last_resort"
        state = "top temperature"
    elif spread_excess > 0.0:
        stage = "last_resort"
        state = "temperature spread"
    else:
        stage = "last_resort"
        state = "nonpositive cooling"
    return {
        "assist_blower_pressure_excess_Pa": assist_blower_pressure_excess,
        "pressure_drop_excess_Pa": pressure_excess,
        "water_loss_excess_kg_day": water_excess,
        "hole_velocity_excess_m_s": hole_velocity_excess,
        "bottom_temp_excess_C": bottom_excess,
        "top_temp_excess_C": top_excess,
        "temp_spread_excess_C": spread_excess,
        "max_bed_temp_C": max(bottom, top),
        "mean_bed_temp_C": 0.5 * (bottom + top),
        "waterLoss_kg_day": water_loss,
        "spot_cooler_load_W": spot_cooler_load,
        "spot_cooler_power_W": spot_cooler_power,
        "assist_blower_power_W": assist_blower_power,
        "combined_cooling_power_W": combined_cooling_power,
        "totalFlow_Lpm": float(pt.get("totalFlow_Lpm", np.nan)),
        "temp_spread_C": spread,
        "minus_cooling_capacity_W": -cooling_capacity,
        "bottom_temp_C": bottom,
        "top_temp_C": top,
        "selectionStage": stage,
        "constraintState": state,
        "isConstraintSafeByCriteria": hard_safe,
        "isFeasibleByCriteria": hard_safe,
    }


def cooling_criteria_sort_key(metrics: dict[str, float | str], config: dict) -> tuple[float, ...]:
    return tuple(float(metrics.get(name, np.inf)) for name in cooling_mode_config(config)["optimization"]["priority_order"])


def rank_cooling_points(payload: dict, config: dict, points: list[dict] | None = None) -> list[dict]:
    source = payload["design_points"] if points is None else points
    ranked: list[dict] = []
    for pt in source:
        metrics = cooling_point_metrics(pt, config)
        enriched = dict(pt)
        enriched["criteriaKey"] = cooling_criteria_sort_key(metrics, config)
        enriched["selectionStage"] = str(metrics["selectionStage"])
        enriched["constraintState"] = str(metrics["constraintState"])
        enriched["isConstraintSafeByCriteria"] = bool(metrics["isConstraintSafeByCriteria"])
        enriched["isFeasibleByCriteria"] = bool(metrics["isFeasibleByCriteria"])
        enriched["tempSpread_C"] = float(metrics["temp_spread_C"])
        enriched["TbottomFullPower_C"] = float(metrics["bottom_temp_C"])
        enriched["TtopFullPower_C"] = float(metrics["top_temp_C"])
        ranked.append(enriched)
    ranked.sort(key=lambda pt: pt["criteriaKey"])
    return ranked


def select_recommended_cooling_point(payload: dict, config: dict, points: list[dict] | None = None) -> dict:
    ranked = rank_cooling_points(payload, config, points)
    hard_safe = [pt for pt in ranked if bool(pt.get("isConstraintSafeByCriteria", False))]
    return hard_safe[0] if hard_safe else {}


def select_best_available_cooling_point(payload: dict, config: dict, points: list[dict] | None = None) -> dict:
    ranked = rank_cooling_points(payload, config, points)
    return ranked[0] if ranked else {}


def cooling_top_candidates(payload: dict, config: dict) -> list[dict]:
    limit = int(cooling_mode_config(config)["optimization"].get("top_candidate_count", 12))
    ranked = rank_cooling_points(payload, config)
    return ranked[: min(limit, len(ranked))]


def cooling_optimization_ranked_points(payload: dict, config: dict) -> list[dict]:
    opt = cooling_mode_config(config)["optimization"]
    if not opt.get("enabled", True):
        return []
    ranked = rank_cooling_points(payload, config)
    return ranked[: int(opt.get("n_top", 10))]


def build_summer_cooling_payload(config: dict) -> dict:
    summary_inputs = config.get("summary_inputs", {})
    model_overrides = normalized_model_overrides(config)
    cooling_cfg = cooling_mode_config(config)
    parallel_cfg = python_parallel_config(config)
    env_inputs = effective_environment_inputs(summary_inputs, model_overrides)
    env_inputs["greenhouseAir_C"] = float(cooling_cfg["ambientAir_C"])
    env_inputs["pressure_Pa"] = float(cooling_cfg.get("pressure_Pa", env_inputs.get("pressure_Pa", 101325.0)))
    design_bed_C = float(cooling_cfg["designBed_C"])
    bin_inputs = effective_bin_inputs(summary_inputs, model_overrides)
    wall_inputs = effective_bin_wall_inputs(summary_inputs, model_overrides)
    bin_model = build_bin_model_from_inputs(bin_inputs, wall_inputs, env_inputs, design_bed_C)
    lumped_loss = lumped_loss_model_from_bin(bin_model, float(env_inputs["greenhouseAir_C"]), design_bed_C)
    required_cooling_W = max(float(env_inputs["greenhouseAir_C"]) - design_bed_C, 0.0) * float(bin_model["UAtotalAmbient_W_K"])
    cooling_flows = cooling_flow_array(cooling_cfg)
    use_parallel = bool(parallel_cfg.get("enabled", False) and parallel_cfg.get("coolingSweep", True) and cooling_flows.size > 1)
    if use_parallel:
        ctx = mp.get_context("spawn")
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=python_process_pool_workers(parallel_cfg),
            mp_context=ctx,
        ) as executor:
            tasks = [(config, float(flow), bin_model, required_cooling_W) for flow in cooling_flows]
            points = list(executor.map(simulate_summer_cooling_point_task, tasks))
    else:
        points = [simulate_summer_cooling_point(config, float(flow), bin_model, required_cooling_W) for flow in cooling_flows]
    payload = {
        "meta": {
            "operatingMode": COOLING_MODE,
            "recommended_definition": (
                "Recommended = the single hard-constraint-safe summer plenum-fed spot-cooler plus assist-blower point chosen by "
                "lexicographic criteria. Inside the hard-safe set, the rule is: minimize the hottest bed-node "
                "temperature, then minimize mean bed temperature, then minimize combined cooling-system electrical power, then minimize airflow."
            ),
        },
        "active_configuration": {"label": "Covered top + vent"},
        "design_points": points,
        "recommended": {},
        "best_available": {},
        "limits": {
            "maxAssistBlowerDeliverablePressure_Pa": float(cooling_cfg["limits"]["maxAssistBlowerDeliverablePressure_Pa"]),
            "maxPressureDrop_Pa": float(cooling_cfg["limits"]["maxPressureDrop_Pa"]),
            "maxWaterLoss_kg_day": float(cooling_cfg["limits"]["maxWaterLoss_kg_day"]),
            "maxHoleVelocity_m_s": float(cooling_cfg["limits"]["maxHoleVelocity_m_s"]),
        },
        "requiredCooling_W": required_cooling_W,
        "requiredHeat_W": -required_cooling_W,
        "lumped_loss": {
            "requiredHeat_W": -required_cooling_W,
            "requiredCooling_W": required_cooling_W,
            "UA_W_K": lumped_loss["UA_W_K"],
            "tau_h": lumped_loss["tau_h"],
            "Bi_lumped": lumped_loss["Bi_lumped"],
            "lumpedStrictlyValid": lumped_loss["lumpedStrictlyValid"],
            "externalH_W_m2K": lumped_loss["externalH_W_m2K"],
        },
        "bin_model": bin_model,
        "configuration_comparisons": [],
        "sweep": {
            "voltage_V": [0.0],
            "totalFlow_Lpm": cooling_flows.tolist(),
            "targetDuty": 1.0,
        },
        "minimumFlowAt100Duty_Lpm": float("nan"),
        "minimumFlowAtTargetDuty_Lpm": float("nan"),
    }
    payload["recommended"] = select_recommended_cooling_point(payload, config)
    payload["best_available"] = select_best_available_cooling_point(payload, config)
    return payload


def python_process_pool_workers(parallel_cfg: dict) -> int | None:
    workers = int(parallel_cfg.get("workers", 0) or 0)
    return None if workers <= 0 else workers


def simulate_summer_cooling_point_task(task: tuple[dict, float, dict, float]) -> dict:
    config, flow_lpm, bin_model, required_cooling_W = task
    return simulate_summer_cooling_point(config, float(flow_lpm), bin_model, float(required_cooling_W))


def compute_wire_length_m(pitch_mm: float, heater: dict, wire: dict) -> float | None:
    coil_span_m = heater.get("coilSpan_m")
    coil_mean_d_mm = heater.get("coilMeanD_mm")
    wire_diameter_mm = wire.get("diameter_mm")
    if coil_span_m is None or coil_mean_d_mm is None:
        return None
    pitch_m = pitch_mm / 1000.0
    if wire_diameter_mm is not None:
        pitch_m = max(pitch_m, 1.05 * wire_diameter_mm / 1000.0)
    turns = max(1.0, coil_span_m / max(pitch_m, 1e-9))
    turn_length_m = float(np.hypot(np.pi * (coil_mean_d_mm / 1000.0), pitch_m))
    return turns * turn_length_m


def bed_fill_volume_m3(summary_inputs: dict, model_overrides: dict | None = None) -> float:
    if model_overrides is None:
        bin_inputs = summary_inputs.get("bin", {})
    else:
        bin_inputs = effective_bin_inputs(summary_inputs, model_overrides)
    required = {"length_m", "width_m", "height_m", "fillFraction"}
    if required <= bin_inputs.keys():
        return (
            float(bin_inputs["length_m"])
            * float(bin_inputs["width_m"])
            * float(bin_inputs["height_m"])
            * float(bin_inputs["fillFraction"])
        )
    return float("nan")


def particle_surface_details(summary_inputs: dict, model_overrides: dict) -> dict:
    evaporation_inputs = model_overrides.get("evaporation", {})
    table3 = dict(evaporation_inputs.get("table3ParticleSurface", {}))
    representatives = dict(table3.get("representativeDiameters_mm", {}))
    windrow = dict(table3.get("windrowFractions_pct", {}))
    vermi = dict(table3.get("vermicompostFractions_pct", {}))
    blend_mode = str(table3.get("blendMode", "mean_of_windrow_and_vermicompost")).strip().lower()
    accessible_fraction = float(
        evaporation_inputs.get(
            "accessibleSurfaceFraction",
            evaporation_inputs.get("wettedAreaFraction", 1.0),
        )
    )

    representative_mm = np.array(
        [float(representatives.get(key, np.nan)) for key in PARTICLE_SURFACE_KEYS],
        dtype=float,
    )
    windrow_pct = np.array([float(windrow.get(key, 0.0)) for key in PARTICLE_SURFACE_KEYS], dtype=float)
    vermi_pct = np.array([float(vermi.get(key, 0.0)) for key in PARTICLE_SURFACE_KEYS], dtype=float)
    if blend_mode == "mean_of_windrow_and_vermicompost":
        blended_pct = 0.5 * (windrow_pct + vermi_pct)
    else:
        blended_pct = windrow_pct.copy()

    weight_sum = float(np.sum(np.clip(blended_pct, 0.0, None)))
    if weight_sum <= 0 or not np.all(np.isfinite(representative_mm)):
        return {
            "basis": "table3_particle_surface",
            "accessible_fraction": max(0.0, accessible_fraction),
            "blended_weights": np.full(len(PARTICLE_SURFACE_KEYS), np.nan),
            "windrow_pct": windrow_pct,
            "vermi_pct": vermi_pct,
            "representative_mm": representative_mm,
            "d32_m": np.nan,
            "area_per_volume_m2_m3": np.nan,
            "fill_volume_m3": bed_fill_volume_m3(summary_inputs, model_overrides),
            "area_m2": np.nan,
        }

    blended_weights = np.clip(blended_pct, 0.0, None) / weight_sum
    representative_m = representative_mm / 1000.0
    inv_d32 = float(np.sum(blended_weights / np.maximum(representative_m, 1e-12)))
    d32_m = 1.0 / max(inv_d32, 1e-12)
    fill_volume_m3 = bed_fill_volume_m3(summary_inputs, model_overrides)
    area_per_volume = 6.0 / max(d32_m, 1e-12)
    area_m2 = max(0.0, accessible_fraction) * fill_volume_m3 * area_per_volume

    return {
        "basis": "table3_particle_surface",
        "accessible_fraction": max(0.0, accessible_fraction),
        "blended_weights": blended_weights,
        "windrow_pct": windrow_pct,
        "vermi_pct": vermi_pct,
        "representative_mm": representative_mm,
        "d32_m": d32_m,
        "area_per_volume_m2_m3": area_per_volume,
        "fill_volume_m3": fill_volume_m3,
        "area_m2": area_m2,
    }


def derived_wetted_area_details(summary_inputs: dict, model_overrides: dict) -> dict:
    evaporation_inputs = model_overrides.get("evaporation", {})
    basis = str(evaporation_inputs.get("wettedAreaBasis", "manual")).strip().lower()
    if basis == "manual" and evaporation_inputs.get("wettedArea_m2") is not None:
        return {
            "basis": basis,
            "area_m2": float(evaporation_inputs["wettedArea_m2"]),
            "fill_volume_m3": bed_fill_volume_m3(summary_inputs, model_overrides),
        }
    if basis == "table3_particle_surface":
        return particle_surface_details(summary_inputs, model_overrides)
    bin_inputs = effective_bin_inputs(summary_inputs, model_overrides)
    if {"length_m", "width_m"} <= bin_inputs.keys():
        wet_fraction = float(evaporation_inputs.get("wettedAreaFraction", 1.0))
        return {
            "basis": basis,
            "wet_fraction": max(0.0, wet_fraction),
            "area_m2": max(0.0, wet_fraction)
            * float(bin_inputs["length_m"])
            * float(bin_inputs["width_m"]),
            "fill_volume_m3": bed_fill_volume_m3(summary_inputs, model_overrides),
        }
    return {
        "basis": basis,
        "area_m2": float(evaporation_inputs.get("wettedArea_m2", np.nan)),
        "fill_volume_m3": bed_fill_volume_m3(summary_inputs, model_overrides),
    }


def derived_wetted_area_m2(summary_inputs: dict, model_overrides: dict) -> float:
    return float(derived_wetted_area_details(summary_inputs, model_overrides).get("area_m2", np.nan))


def effective_bin_inputs(summary_inputs: dict, model_overrides: dict) -> dict:
    bin_inputs = dict(summary_inputs.get("bin", {}))
    bin_inputs.update(model_overrides.get("bin", {}))
    return bin_inputs


def effective_environment_inputs(summary_inputs: dict, model_overrides: dict) -> dict:
    env_inputs = {
        "greenhouseAir_C": float(summary_inputs.get("greenhouseAir_C", -15.0)),
        "pressure_Pa": float(summary_inputs.get("pressure_Pa", 101325.0)),
        "referencePlateLength_m": float(
            summary_inputs.get(
                "referencePlateLength_m",
                summary_inputs.get("bin", {}).get("height_m", 0.61),
            )
        ),
    }
    env_inputs.update(model_overrides.get("environment", {}))
    return env_inputs


def effective_bin_wall_inputs(summary_inputs: dict, model_overrides: dict) -> dict:
    wall_inputs = dict(summary_inputs.get("bin_wall", {}))
    wall_override = dict(model_overrides.get("binWall", {}))
    if "sheetThickness_m" in wall_override:
        wall_override["sheetThickness_mm"] = 1000.0 * float(wall_override["sheetThickness_m"])
    if "insulationThickness_m" in wall_override:
        wall_override["insulationThickness_mm"] = 1000.0 * float(wall_override["insulationThickness_m"])
    wall_inputs.update(wall_override)
    return wall_inputs


def effective_heater_inputs(summary_inputs: dict, model_overrides: dict) -> dict:
    heater_inputs = dict(summary_inputs.get("heater_tube", {}))
    heater_override = dict(model_overrides.get("heaterTube", {}))
    if "ID_m" in heater_override:
        heater_override["ID_mm"] = 1000.0 * float(heater_override["ID_m"])
    if "OD_m" in heater_override:
        heater_override["OD_mm"] = 1000.0 * float(heater_override["OD_m"])
    if "insulationThickness_m" in heater_override:
        heater_override["insulationThickness_mm"] = 1000.0 * float(heater_override["insulationThickness_m"])
    if "coilMeanD_m" in heater_override:
        heater_override["coilMeanD_mm"] = 1000.0 * float(heater_override["coilMeanD_m"])
    if "coilSpan_m" in heater_override:
        heater_override["coilSpan_m"] = float(heater_override["coilSpan_m"])
    heater_inputs.update(heater_override)
    return heater_inputs


def effective_aeration_inputs(summary_inputs: dict, model_overrides: dict) -> dict:
    aer_inputs = dict(summary_inputs.get("aeration", {}))
    aer_overrides = dict(model_overrides.get("aeration", {}))
    if "ID_m" in aer_overrides:
        aer_overrides["ID_mm"] = 1000.0 * float(aer_overrides["ID_m"])
    if "OD_m" in aer_overrides:
        aer_overrides["OD_mm"] = 1000.0 * float(aer_overrides["OD_m"])
    if "headerID_m" in aer_overrides:
        aer_overrides["headerID_mm"] = 1000.0 * float(aer_overrides["headerID_m"])
    if "splitterInletID_m" in aer_overrides:
        aer_overrides["splitterInletID_mm"] = 1000.0 * float(aer_overrides["splitterInletID_m"])
    if "branchConnectorID_m" in aer_overrides:
        aer_overrides["branchConnectorID_mm"] = 1000.0 * float(aer_overrides["branchConnectorID_m"])
    aer_inputs.update(aer_overrides)
    if "nParallelTubes" in aer_inputs:
        aer_inputs["splitterOutletCount"] = max(1, int(round(float(aer_inputs["nParallelTubes"]))))
    return aer_inputs


def effective_wire_inputs(summary_inputs: dict, model_overrides: dict) -> dict:
    wire_inputs = dict(summary_inputs.get("wire", {}))
    wire_overrides = dict(model_overrides.get("wire", {}))
    if "diameter_m" in wire_overrides:
        wire_overrides["diameter_mm"] = 1000.0 * float(wire_overrides["diameter_m"])
    wire_inputs.update(wire_overrides)
    wire_inputs.setdefault("rho20_Ohm_m", 1.09e-6)
    wire_inputs.setdefault("tempCoeff_1_K", 1.7e-4)
    return wire_inputs


def build_heating_wire_geometry(point: dict, heater_inputs: dict, wire_inputs: dict) -> dict:
    wire_diameter_m = float(wire_inputs.get("diameter_mm", np.nan)) / 1000.0
    coil_span_m = float(heater_inputs.get("coilSpan_m", heater_inputs.get("length_m", np.nan)))
    coil_mean_d_m = float(heater_inputs.get("coilMeanD_mm", np.nan)) / 1000.0
    heater_id_m = float(heater_inputs.get("ID_mm", np.nan)) / 1000.0
    pitch_m = float(point.get("pitch_mm", np.nan)) / 1000.0
    if np.isfinite(wire_diameter_m):
        pitch_m = max(pitch_m, 1.05 * wire_diameter_m)

    is_valid = all(
        np.isfinite(val) and val > 0.0
        for val in (wire_diameter_m, coil_span_m, coil_mean_d_m, heater_id_m, pitch_m)
    )
    fit_clearance_m = float("nan")
    if is_valid:
        fit_clearance_m = heater_id_m - (coil_mean_d_m + wire_diameter_m)
        is_valid = fit_clearance_m >= 0.0

    wire_length_m = float(point.get("wireLength_m", np.nan))
    if not np.isfinite(wire_length_m) or wire_length_m <= 0.0:
        turn_length_m = float(np.hypot(np.pi * coil_mean_d_m, pitch_m)) if is_valid else float("nan")
        turns = max(coil_span_m / max(pitch_m, 1e-12), 1.0) if is_valid else float("nan")
        wire_length_m = turns * turn_length_m if is_valid else float("nan")
    else:
        turn_length_m = float(np.hypot(np.pi * coil_mean_d_m, pitch_m)) if is_valid else float("nan")
        turns = wire_length_m / max(turn_length_m, 1e-12) if is_valid else float("nan")

    wire_cross_section_m2 = np.pi * wire_diameter_m**2 / 4.0 if np.isfinite(wire_diameter_m) else float("nan")
    wire_surface_area_m2 = np.pi * wire_diameter_m * wire_length_m if np.isfinite(wire_length_m) else float("nan")
    wire_surface_area_per_axial_length_m = wire_surface_area_m2 / max(coil_span_m, 1e-12)

    return {
        "isValid": bool(is_valid),
        "fitClearance_m": float(fit_clearance_m),
        "outerDiameter_m": float(coil_mean_d_m + wire_diameter_m) if np.isfinite(coil_mean_d_m) and np.isfinite(wire_diameter_m) else float("nan"),
        "turns": float(turns) if np.isfinite(turns) else float("nan"),
        "pitch_m": float(pitch_m) if np.isfinite(pitch_m) else float("nan"),
        "length_m": float(wire_length_m) if np.isfinite(wire_length_m) else float("nan"),
        "crossSection_m2": float(wire_cross_section_m2) if np.isfinite(wire_cross_section_m2) else float("nan"),
        "surfaceArea_m2": float(wire_surface_area_m2) if np.isfinite(wire_surface_area_m2) else float("nan"),
        "surfaceAreaPerAxialLength_m": float(wire_surface_area_per_axial_length_m) if np.isfinite(wire_surface_area_per_axial_length_m) else float("nan"),
        "turnLength_m": float(turn_length_m) if np.isfinite(turn_length_m) else float("nan"),
        "coilMeanD_m": float(coil_mean_d_m) if np.isfinite(coil_mean_d_m) else float("nan"),
        "diameter_m": float(wire_diameter_m) if np.isfinite(wire_diameter_m) else float("nan"),
    }


def normalize_string_tube_counts(raw_counts, n_tubes: int) -> list[int] | None:
    if raw_counts is None:
        return None
    counts: list[int] = []
    try:
        if isinstance(raw_counts, str):
            cleaned = raw_counts.strip().strip("[]")
            if not cleaned:
                return None
            tokens = cleaned.replace(",", " ").split()
            counts = [int(round(float(token))) for token in tokens]
        else:
            counts = [int(round(float(value))) for value in raw_counts]
    except (TypeError, ValueError):
        return None
    counts = [count for count in counts if count > 0]
    if not counts or sum(counts) != max(1, int(n_tubes)):
        return None
    return counts


def derive_heating_string_tube_counts(n_tubes: int, electrical_cfg: dict | None = None) -> list[int]:
    n_tubes = max(1, int(round(n_tubes)))
    electrical_cfg = dict(electrical_cfg or {})
    manual = normalize_string_tube_counts(electrical_cfg.get("manualStringTubeCounts"), n_tubes)
    if manual is not None:
        return manual
    mode = str(electrical_cfg.get("topologyMode", "paired_series_strings")).strip().lower()
    if mode == "all_parallel":
        return [1] * n_tubes
    if n_tubes == 1:
        return [1]
    if n_tubes % 2 == 0:
        return [2] * (n_tubes // 2)
    return [3] + [2] * ((n_tubes - 3) // 2)


def heating_string_tube_counts(model_overrides: dict, point: dict, n_tubes: int) -> list[int]:
    exported = normalize_string_tube_counts(point.get("stringTubeCounts"), n_tubes)
    if exported is not None:
        return exported
    electrical_cfg = model_overrides.get("electrical", {})
    return derive_heating_string_tube_counts(n_tubes, electrical_cfg)


def electrical_topology_label(string_counts: list[int]) -> str:
    if not string_counts:
        return "unspecified"
    parts: list[str] = []
    for n_series in sorted(set(string_counts)):
        n_strings = sum(1 for value in string_counts if value == n_series)
        parts.append(f"{n_strings}x{n_series}S")
    return " + ".join(parts)


def format_string_tube_counts(string_counts: list[int]) -> str:
    return "[" + ", ".join(str(int(count)) for count in string_counts) + "]"


def air_props(T_C: float, pressure_Pa: float = 101325.0) -> dict:
    T_K = T_C + 273.15
    rho = pressure_Pa / max(287.058 * T_K, 1e-12)
    mu0 = 1.716e-5
    T0 = 273.15
    suth = 111.0
    mu = mu0 * (T_K / T0) ** 1.5 * (T0 + suth) / max(T_K + suth, 1e-12)
    cp = 1007.0
    pr = 0.71
    k = mu * cp / pr
    return {"rho_kg_m3": rho, "mu_Pa_s": mu, "cp_J_kgK": cp, "Pr": pr, "k_W_mK": k}


def saturation_pressure_water_Pa(T_K: float) -> float:
    Tc_K = 647.096
    pc_Pa = 22.064e6
    tau = 1.0 - T_K / Tc_K
    coeffs = (-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502)
    exponents = (1.0, 1.5, 3.0, 3.5, 4.0, 7.5)
    series = sum(a * tau**b for a, b in zip(coeffs, exponents))
    return pc_Pa * np.exp((Tc_K / max(T_K, 1e-12)) * series)


def saturation_pressure_slope_water_Pa_K(T_K: float) -> float:
    Tc_K = 647.096
    pc_Pa = 22.064e6
    tau = 1.0 - T_K / Tc_K
    coeffs = (-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502)
    exponents = (1.0, 1.5, 3.0, 3.5, 4.0, 7.5)
    sum_series = sum(a * tau**b for a, b in zip(coeffs, exponents))
    sum_deriv = sum(a * b * tau ** (b - 1.0) for a, b in zip(coeffs, exponents))
    lnp = (Tc_K / max(T_K, 1e-12)) * sum_series
    psat = pc_Pa * np.exp(lnp)
    dln_dT = -(Tc_K / max(T_K, 1e-12) ** 2) * sum_series - sum_deriv / max(T_K, 1e-12)
    return psat * dln_dT


def saturation_liquid_density_water_kg_m3(T_K: float) -> float:
    Tc_K = 647.096
    rhoc_kg_m3 = 322.0
    tau = 1.0 - T_K / Tc_K
    coeffs = (1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -6.74694450e5)
    exponents = (1.0 / 3.0, 2.0 / 3.0, 5.0 / 3.0, 16.0 / 3.0, 43.0 / 3.0, 110.0 / 3.0)
    series = 1.0 + sum(a * tau**b for a, b in zip(coeffs, exponents))
    return rhoc_kg_m3 * series


def saturation_vapor_density_water_kg_m3(T_K: float) -> float:
    Tc_K = 647.096
    rhoc_kg_m3 = 322.0
    tau = 1.0 - T_K / Tc_K
    coeffs = (-2.03150240, -2.68302940, -5.38626492, -17.2991605, -44.7586581, -63.9201063)
    exponents = (2.0 / 6.0, 4.0 / 6.0, 8.0 / 6.0, 18.0 / 6.0, 37.0 / 6.0, 71.0 / 6.0)
    series = sum(a * tau**b for a, b in zip(coeffs, exponents))
    return rhoc_kg_m3 * np.exp(series)


def water_vapor_density_at_saturation_kg_m3(T_K: float) -> float:
    Rv_J_kgK = 461.526
    return saturation_pressure_water_Pa(T_K) / max(Rv_J_kgK * T_K, 1e-12)


def latent_heat_vaporization_water_J_kg(T_K: float) -> float:
    rho_l = saturation_liquid_density_water_kg_m3(T_K)
    rho_v = saturation_vapor_density_water_kg_m3(T_K)
    dpdT = saturation_pressure_slope_water_Pa_K(T_K)
    return T_K * dpdT * (1.0 / max(rho_v, 1e-12) - 1.0 / max(rho_l, 1e-12))


def humidity_ratio_from_rh(T_C: float, relative_humidity: float, pressure_Pa: float) -> float:
    rh = min(max(relative_humidity, 0.0), 1.0)
    p_sat = saturation_pressure_water_Pa(T_C + 273.15)
    p_v = min(rh * p_sat, 0.999 * pressure_Pa)
    return 0.62198 * p_v / max(pressure_Pa - p_v, 1e-12)


def humidity_ratio_saturated(T_C: float, pressure_Pa: float) -> float:
    p_sat = min(saturation_pressure_water_Pa(T_C + 273.15), 0.999 * pressure_Pa)
    return 0.62198 * p_sat / max(pressure_Pa - p_sat, 1e-12)


def relative_humidity_from_humidity_ratio(T_C: float, humidity_ratio: float, pressure_Pa: float) -> float:
    w = max(humidity_ratio, 0.0)
    p_v = pressure_Pa * w / max(0.62198 + w, 1e-12)
    p_sat = saturation_pressure_water_Pa(T_C + 273.15)
    return min(max(p_v / max(p_sat, 1e-12), 0.0), 1.0)


def dew_point_from_rh(T_C: float, relative_humidity: float, pressure_Pa: float) -> float:
    rh = min(max(relative_humidity, 0.0), 1.0)
    p_target = rh * saturation_pressure_water_Pa(T_C + 273.15)
    lo_C = -40.0
    hi_C = T_C
    for _ in range(60):
        mid_C = 0.5 * (lo_C + hi_C)
        p_mid = saturation_pressure_water_Pa(mid_C + 273.15)
        if p_mid < p_target:
            lo_C = mid_C
        else:
            hi_C = mid_C
    return 0.5 * (lo_C + hi_C)


def natural_convection_vertical_plate_h(Tsurf_C: float, Tinf_C: float, L_m: float, pressure_Pa: float) -> float:
    Tfilm_K = 0.5 * (Tsurf_C + Tinf_C) + 273.15
    props = air_props(0.5 * (Tsurf_C + Tinf_C), pressure_Pa)
    beta_1_K = 1.0 / max(Tfilm_K, 1e-12)
    nu_m2_s = props["mu_Pa_s"] / max(props["rho_kg_m3"], 1e-12)
    alpha_m2_s = props["k_W_mK"] / max(props["rho_kg_m3"] * props["cp_J_kgK"], 1e-12)
    deltaT_K = max(abs(Tsurf_C - Tinf_C), 1e-6)
    Ra = 9.81 * beta_1_K * deltaT_K * max(L_m, 1e-9) ** 3 / max(nu_m2_s * alpha_m2_s, 1e-12)
    Nu = (0.825 + 0.387 * Ra ** (1.0 / 6.0) / (1.0 + (0.492 / max(props["Pr"], 1e-12)) ** (9.0 / 16.0)) ** (8.0 / 27.0)) ** 2
    return Nu * props["k_W_mK"] / max(L_m, 1e-9)


def natural_convection_horizontal_plate_h(Tsurf_C: float, Tinf_C: float, L_m: float, pressure_Pa: float) -> float:
    Tfilm_K = 0.5 * (Tsurf_C + Tinf_C) + 273.15
    props = air_props(0.5 * (Tsurf_C + Tinf_C), pressure_Pa)
    beta_1_K = 1.0 / max(Tfilm_K, 1e-12)
    nu_m2_s = props["mu_Pa_s"] / max(props["rho_kg_m3"], 1e-12)
    alpha_m2_s = props["k_W_mK"] / max(props["rho_kg_m3"] * props["cp_J_kgK"], 1e-12)
    deltaT_K = max(abs(Tsurf_C - Tinf_C), 1e-6)
    Ra = 9.81 * beta_1_K * deltaT_K * max(L_m, 1e-9) ** 3 / max(nu_m2_s * alpha_m2_s, 1e-12)
    if Tsurf_C >= Tinf_C:
        Nu = 0.54 * Ra ** 0.25 if Ra <= 1.0e7 else 0.15 * Ra ** (1.0 / 3.0)
    else:
        Nu = 0.27 * Ra ** 0.25
    return Nu * props["k_W_mK"] / max(L_m, 1e-9)


def overall_wall_u(wall_inputs: dict, h_outside_W_m2K: float) -> float:
    rin = 1.0 / max(float(wall_inputs.get("internalH_W_m2K", 6.0)), 1e-12)
    rsheet = float(wall_inputs.get("sheetThickness_mm", 1.5)) / 1000.0 / max(float(wall_inputs.get("sheetK_W_mK", 50.0)), 1e-12)
    rins = float(wall_inputs.get("insulationThickness_mm", 50.8)) / 1000.0 / max(float(wall_inputs.get("insulationK_W_mK", 0.040)), 1e-12)
    rout = 1.0 / max(h_outside_W_m2K, 1e-12)
    return 1.0 / max(rin + rsheet + rins + rout, 1e-12)


def churchill_bernstein_nu(Re: float, Pr: float) -> float:
    Re = max(Re, 1e-12)
    Pr = max(Pr, 1e-12)
    return 0.3 + (
        (0.62 * Re**0.5 * Pr ** (1.0 / 3.0))
        / (1.0 + (0.4 / Pr) ** (2.0 / 3.0)) ** 0.25
        * (1.0 + (Re / 282000.0) ** (5.0 / 8.0)) ** (4.0 / 5.0)
    )


def churchill_friction_factor(Re: float) -> float:
    Re = max(Re, 1e-12)
    if Re < 2300.0:
        return 64.0 / Re
    A = (2.457 * np.log(1.0 / ((7.0 / Re) ** 0.9))) ** 16
    B = (37530.0 / Re) ** 16
    return 8.0 * ((8.0 / Re) ** 12 + 1.0 / (A + B) ** 1.5) ** (1.0 / 12.0)


def internal_tube_nu(Re: float, Pr: float, D_m: float, L_m: float) -> float:
    Re = max(Re, 1e-12)
    Pr = max(Pr, 1e-12)
    if Re < 2300.0:
        gz = Re * Pr * max(D_m, 1e-12) / max(L_m, 1e-12)
        return 3.66 + (0.0668 * gz) / (1.0 + 0.04 * gz ** (2.0 / 3.0))
    f = churchill_friction_factor(Re)
    nu = (f / 8.0) * (Re - 1000.0) * Pr / (1.0 + 12.7 * (f / 8.0) ** 0.5 * (Pr ** (2.0 / 3.0) - 1.0))
    return max(nu, 3.66)


def heater_tube_overall_u(heater_inputs: dict, h_inside_wall_W_m2K: float) -> float:
    rin = 1.0 / max(h_inside_wall_W_m2K, 1e-12)
    rwall = (
        (float(heater_inputs.get("OD_mm", 19.0)) - float(heater_inputs.get("ID_mm", 15.0))) / 1000.0
    ) / max(2.0 * float(heater_inputs.get("wallK_W_mK", 16.0)), 1e-12)
    rins = (float(heater_inputs.get("insulationThickness_mm", 10.0)) / 1000.0) / max(
        float(heater_inputs.get("insulationK_W_mK", 0.040)), 1e-12
    )
    rout = 1.0 / max(float(heater_inputs.get("externalH_W_m2K", 10.0)), 1e-12)
    return 1.0 / max(rin + rwall + rins + rout, 1e-12)


def solve_heater_tube_python(
    heater_inputs: dict,
    wire_inputs: dict,
    calibration: dict,
    pressure_Pa: float,
    wire: dict,
    power_W: float,
    flow_Lpm: float,
    air_in_K: float,
    ambient_C: float,
) -> dict:
    n_segments = max(1, int(round(float(heater_inputs.get("segments", 50)))))
    length_m = max(float(heater_inputs.get("length_m", 0.30)), 1e-12)
    dx_m = length_m / n_segments
    id_m = max(float(heater_inputs.get("ID_mm", 15.0)) / 1000.0, 1e-12)
    od_m = max(float(heater_inputs.get("OD_mm", 19.0)) / 1000.0, id_m)
    atube_m2 = np.pi * id_m**2 / 4.0
    pout_m = np.pi * od_m
    tair_K = float(air_in_K)

    props_in = air_props(air_in_K - 273.15, pressure_Pa)
    vol_flow_m3_s = flow_Lpm / 60000.0
    mdot_kg_s = vol_flow_m3_s * props_in["rho_kg_m3"]

    wire_profile_C = np.zeros(n_segments, dtype=float)
    air_profile_C = np.zeros(n_segments + 1, dtype=float)
    air_profile_C[0] = air_in_K - 273.15
    hwire_profile = np.zeros(n_segments, dtype=float)
    dp_friction_pa = 0.0

    coil_mean_d_m = float(wire.get("coilMeanD_m", np.nan))
    wire_diameter_m = float(wire.get("diameter_m", np.nan))
    pitch_m = float(wire.get("pitch_m", np.nan))

    for idx in range(n_segments):
        props = air_props(tair_K - 273.15, pressure_Pa)
        u_bulk_m_s = mdot_kg_s / max(props["rho_kg_m3"] * atube_m2, 1e-12)
        helix_angle_factor = (np.pi * coil_mean_d_m) / max(np.hypot(np.pi * coil_mean_d_m, pitch_m), 1e-12)
        u_normal_m_s = max(1e-6, u_bulk_m_s * helix_angle_factor)
        re_wire = props["rho_kg_m3"] * u_normal_m_s * wire_diameter_m / max(props["mu_Pa_s"], 1e-12)
        nu_wire = churchill_bernstein_nu(re_wire, props["Pr"])
        h_wire = float(calibration.get("hWireMultiplier", 1.0)) * nu_wire * props["k_W_mK"] / max(wire_diameter_m, 1e-12)

        re_tube = props["rho_kg_m3"] * u_bulk_m_s * id_m / max(props["mu_Pa_s"], 1e-12)
        nu_tube = internal_tube_nu(re_tube, props["Pr"], id_m, length_m)
        h_inside_wall = float(calibration.get("hWallMultiplier", 1.0)) * nu_tube * props["k_W_mK"] / id_m
        u_loss = heater_tube_overall_u(heater_inputs, h_inside_wall)

        q_loss_W = u_loss * pout_m * dx_m * max((tair_K - 273.15) - ambient_C, 0.0)
        q_gen_W = power_W * dx_m / max(length_m, 1e-12)
        q_to_air_W = q_gen_W - q_loss_W
        tair_K = tair_K + q_to_air_W / max(mdot_kg_s * props["cp_J_kgK"], 1e-12)

        delta_twire_K = q_gen_W / max(h_wire * float(wire.get("surfaceAreaPerAxialLength_m", np.nan)) * dx_m, 1e-12)
        wire_profile_C[idx] = (tair_K - 273.15) + delta_twire_K
        air_profile_C[idx + 1] = tair_K - 273.15
        hwire_profile[idx] = h_wire

        f = churchill_friction_factor(re_tube)
        dp_friction_pa += (
            f * dx_m / id_m + float(heater_inputs.get("extraMinorK", 1.0)) / n_segments
        ) * 0.5 * props["rho_kg_m3"] * u_bulk_m_s**2

    return {
        "airInletTemp_C": float(air_in_K - 273.15),
        "airOutletTemp_C": float(tair_K - 273.15),
        "airProfile_C": air_profile_C,
        "wireProfile_C": wire_profile_C,
        "wireMaxTemp_C": float(np.max(wire_profile_C)),
        "wireMeanTemp_C": float(np.mean(wire_profile_C)),
        "hWireMean_W_m2K": float(np.mean(hwire_profile)),
        "deltaP_Pa": float(dp_friction_pa),
    }


def solve_heating_aeration_tube_python(
    summary_inputs: dict,
    model_overrides: dict,
    env_inputs: dict,
    aer_inputs: dict,
    calibration: dict,
    tair_in_C: float,
    tbed_C: float,
    mdot_in_kg_s: float,
) -> dict:
    n_segments = max(1, int(round(float(aer_inputs.get("segments", 40)))))
    length_m = max(float(aer_inputs.get("length_m", 1.22)), 1e-12)
    dx_m = length_m / n_segments
    id_m = max(float(aer_inputs.get("ID_mm", 40.0)) / 1000.0, 1e-12)
    od_m = max(float(aer_inputs.get("OD_mm", 44.0)) / 1000.0, id_m)
    atube_m2 = np.pi * id_m**2 / 4.0
    pout_m = np.pi * od_m
    wall_k = float(aer_inputs.get("wallK_W_mK", 16.0))
    pressure_pa = float(env_inputs.get("pressure_Pa", 101325.0))
    release_fraction = float(aer_inputs.get("releaseFraction", 0.95))

    remaining_mdot_kg_s = mdot_in_kg_s
    released_total_kg_s = mdot_in_kg_s * release_fraction
    dm_release_kg_s = released_total_kg_s / n_segments
    tair_C = float(tair_in_C)
    q_to_bed_W = 0.0
    latent_evap_W = 0.0
    water_loss_kg_s = 0.0
    dp_pa = 0.0

    air_profile_C = np.zeros(n_segments + 1, dtype=float)
    air_profile_C[0] = tair_C

    evap_inputs = model_overrides.get("evaporation", {})
    n_tubes = max(1, int(round(float(aer_inputs.get("nParallelTubes", 1)))))
    hole_diameter_m = max(1e-12, float(model_overrides.get("perforation", {}).get("holeDiameter_m", 0.005)))
    holes_per_tube = max(1, int(round(float(model_overrides.get("perforation", {}).get("holesPerTube", 60)))))
    props_in = air_props(tair_in_C, pressure_pa)
    vol_release_m3_s = release_fraction * mdot_in_kg_s / max(props_in["rho_kg_m3"], 1e-12)
    flow_per_hole_lpm = vol_release_m3_s * 60000.0 / holes_per_tube
    hole_velocity_m_s = vol_release_m3_s / max(holes_per_tube * np.pi * hole_diameter_m**2 / 4.0, 1e-12)
    seg_area_m2 = derived_wetted_area_m2(summary_inputs, model_overrides) / max(n_tubes * n_segments, 1)
    inlet_relative_humidity = float(evap_inputs.get("relativeHumidity", 0.80))
    surface_temp_C = float(evap_inputs.get("surfaceTemp_C", tbed_C))
    lewis_factor = float(evap_inputs.get("lewisFactor", 1.0))
    evaporation_multiplier = float(calibration.get("evaporationHMultiplier", 1.0))

    for idx in range(n_segments):
        props = air_props(tair_C, pressure_pa)
        u_bulk_m_s = remaining_mdot_kg_s / max(props["rho_kg_m3"] * atube_m2, 1e-12)
        re_tube = props["rho_kg_m3"] * u_bulk_m_s * id_m / max(props["mu_Pa_s"], 1e-12)
        nu_tube = internal_tube_nu(re_tube, props["Pr"], id_m, length_m)
        h_inside = float(calibration.get("hWallMultiplier", 1.0)) * nu_tube * props["k_W_mK"] / id_m
        h_outside = float(calibration.get("bedHTMultiplier", 1.0)) * float(aer_inputs.get("bedH_W_m2K", 8.0))
        u_tube = 1.0 / max(
            1.0 / max(h_inside, 1e-12)
            + (od_m - id_m) / max(2.0 * wall_k, 1e-12)
            + 1.0 / max(h_outside, 1e-12),
            1e-12,
        )

        q_wall_W = u_tube * pout_m * dx_m * (tair_C - tbed_C)
        dm_kg_s = min(dm_release_kg_s, remaining_mdot_kg_s)
        q_jet_W = dm_kg_s * props["cp_J_kgK"] * (tair_C - tbed_C)
        q_avail_W = max(q_jet_W, 0.0)
        m_evap_kg_s, q_evap_W = summer_evaporation_loss_segment(
            dm_kg_s,
            tair_C,
            tbed_C,
            q_avail_W,
            hole_velocity_m_s,
            seg_area_m2,
            pressure_pa,
            inlet_relative_humidity,
            surface_temp_C,
            lewis_factor,
            evaporation_multiplier,
            hole_diameter_m,
            float(env_inputs.get("greenhouseAir_C", tbed_C)),
        )
        q_to_bed_W += q_wall_W + q_jet_W - q_evap_W
        latent_evap_W += q_evap_W
        water_loss_kg_s += m_evap_kg_s

        remaining_after_release_kg_s = max(remaining_mdot_kg_s - dm_kg_s, 1e-9)
        tair_C = tair_C - q_wall_W / max(remaining_after_release_kg_s * props["cp_J_kgK"], 1e-12)
        remaining_mdot_kg_s = remaining_after_release_kg_s
        air_profile_C[idx + 1] = tair_C

        f = churchill_friction_factor(re_tube)
        dp_pa += (
            f * dx_m / id_m + float(aer_inputs.get("endMinorK", 1.0)) / n_segments
        ) * 0.5 * props["rho_kg_m3"] * u_bulk_m_s**2

    return {
        "QtoBed_W": float(q_to_bed_W),
        "airOutlet_C": float(tair_C),
        "airProfile_C": air_profile_C,
        "deltaP_Pa": float(dp_pa),
        "flowPerHole_Lpm": float(flow_per_hole_lpm),
        "holeVelocity_m_s": float(hole_velocity_m_s),
        "waterLoss_kg_day": float(water_loss_kg_s * 86400.0),
        "latentEvap_W": float(latent_evap_W),
    }


def simulate_heating_reference_point(config: dict, point: dict, ambient_air_C: float, bed_temp_C: float) -> dict:
    summary_inputs = config.get("summary_inputs", {})
    model_overrides = normalized_model_overrides(config)
    env_inputs = effective_environment_inputs(summary_inputs, model_overrides)
    env_inputs["greenhouseAir_C"] = float(ambient_air_C)
    heater_inputs = effective_heater_inputs(summary_inputs, model_overrides)
    aer_inputs = effective_aeration_inputs(summary_inputs, model_overrides)
    wire_inputs = effective_wire_inputs(summary_inputs, model_overrides)
    calibration = model_overrides.get("calibration", {})
    pressure_pa = float(env_inputs.get("pressure_Pa", 101325.0))

    n_tubes = max(1, int(round(float(point.get("nParallelTubes", aer_inputs.get("nParallelTubes", 1))))))
    string_counts = heating_string_tube_counts(model_overrides, point, n_tubes)
    topology_label = str(point.get("electricalTopologyLabel") or electrical_topology_label(string_counts))
    n_tubes = sum(string_counts)
    total_flow_lpm = float(point.get("totalFlow_Lpm", 0.0))
    branch_flow_lpm = total_flow_lpm / n_tubes
    voltage_V = float(point.get("voltage_V", 0.0))
    air_in_K = ambient_air_C + 273.15
    wire = build_heating_wire_geometry(point, heater_inputs, wire_inputs)
    if not wire.get("isValid", False):
        return {
            "totalPower_W": float("inf"),
            "totalCurrent_A": float("inf"),
            "QtoBed_W": 0.0,
            "deltaP_Pa": float("inf"),
            "airOutlet_C": float(ambient_air_C),
            "wireMax_C": float("inf"),
            "waterLoss_kg_day": float("inf"),
            "latentEvap_W": float("nan"),
            "holeVelocity_m_s": float("inf"),
            "electricalTopologyLabel": topology_label,
            "stringTubeCounts": string_counts,
        }

    avg_wire_C = ambient_air_C + 20.0
    tube_resistance_ohm = float("nan")
    group_states: list[dict] = []
    unique_series_counts = sorted(set(string_counts))

    for _ in range(8):
        wire_area_m2 = max(float(wire.get("crossSection_m2", np.nan)), 1e-12)
        rho_ohm_m = float(wire_inputs.get("rho20_Ohm_m", 1.09e-6)) * (
            1.0 + float(wire_inputs.get("tempCoeff_1_K", 1.7e-4)) * (avg_wire_C - 20.0)
        )
        tube_resistance_ohm = rho_ohm_m * float(wire.get("length_m", np.nan)) / wire_area_m2
        group_states = []
        weighted_wire_mean_C = 0.0
        total_tube_weight = 0
        for series_count in unique_series_counts:
            n_strings_group = sum(1 for count in string_counts if count == series_count)
            n_tubes_group = n_strings_group * series_count
            tube_current_A = voltage_V / max(series_count * tube_resistance_ohm, 1e-12)
            tube_power_W = tube_current_A**2 * tube_resistance_ohm
            heater_state = solve_heater_tube_python(
                heater_inputs,
                wire_inputs,
                calibration,
                pressure_pa,
                wire,
                tube_power_W,
                branch_flow_lpm,
                air_in_K,
                ambient_air_C,
            )
            group_states.append(
                {
                    "series_count": series_count,
                    "n_strings": n_strings_group,
                    "n_tubes": n_tubes_group,
                    "tube_current_A": tube_current_A,
                    "tube_power_W": tube_power_W,
                    "heater_state": heater_state,
                }
            )
            weighted_wire_mean_C += n_tubes_group * float(np.mean(heater_state["wireProfile_C"]))
            total_tube_weight += n_tubes_group
        new_avg_wire_C = weighted_wire_mean_C / max(total_tube_weight, 1)
        if abs(new_avg_wire_C - avg_wire_C) < 0.05:
            avg_wire_C = new_avg_wire_C
            break
        avg_wire_C = 0.65 * avg_wire_C + 0.35 * new_avg_wire_C

    props_at_inlet = air_props(ambient_air_C, pressure_pa)
    vol_flow_branch_m3_s = branch_flow_lpm / 60000.0
    mdot_branch_kg_s = vol_flow_branch_m3_s * props_at_inlet["rho_kg_m3"]
    distribution = branch_distribution_losses(
        total_flow_lpm / 60000.0,
        branch_flow_lpm / 60000.0,
        ambient_air_C,
        pressure_pa,
        n_tubes,
        aer_inputs,
        max(float(heater_inputs.get("ID_mm", 15.0)) / 1000.0, 1e-12),
    )
    total_current_A = 0.0
    total_power_W = 0.0
    q_to_bed_W = 0.0
    water_loss_kg_day = 0.0
    latent_evap_W = 0.0
    max_branch_delta_p = -np.inf
    max_air_outlet_C = -np.inf
    max_wire_C = -np.inf
    max_hole_velocity_m_s = -np.inf
    max_flow_per_hole_lpm = -np.inf
    weighted_hwire = 0.0
    weighted_wire_mean = 0.0
    weighted_branch_q = 0.0
    tube_current_min_A = float("inf")
    tube_current_max_A = -float("inf")
    tube_current_weighted_A = 0.0
    tube_power_min_W = float("inf")
    tube_power_max_W = -float("inf")
    tube_power_weighted_W = 0.0
    for group_state in group_states:
        heater_state = group_state["heater_state"]
        branch_aeration = solve_heating_aeration_tube_python(
            summary_inputs,
            model_overrides,
            env_inputs,
            aer_inputs,
            calibration,
            float(heater_state["airOutletTemp_C"]),
            bed_temp_C,
            mdot_branch_kg_s,
        )
        dp_expansion_pa = sudden_expansion_loss_Pa(
            mdot_branch_kg_s,
            float(heater_state["airOutletTemp_C"]),
            pressure_pa,
            max(float(heater_inputs.get("ID_mm", 15.0)) / 1000.0, 1e-12),
            max(float(aer_inputs.get("ID_mm", 40.0)) / 1000.0, 1e-12),
        )
        branch_delta_p = float(heater_state["deltaP_Pa"]) + dp_expansion_pa + float(branch_aeration["deltaP_Pa"])
        n_tubes_group = int(group_state["n_tubes"])
        n_strings_group = int(group_state["n_strings"])
        tube_current_A = float(group_state["tube_current_A"])
        tube_power_W = float(group_state["tube_power_W"])

        total_current_A += n_strings_group * tube_current_A
        total_power_W += n_tubes_group * tube_power_W
        q_to_bed_W += n_tubes_group * float(branch_aeration["QtoBed_W"])
        water_loss_kg_day += n_tubes_group * float(branch_aeration["waterLoss_kg_day"])
        latent_evap_W += n_tubes_group * float(branch_aeration["latentEvap_W"])
        max_branch_delta_p = max(max_branch_delta_p, branch_delta_p)
        max_air_outlet_C = max(max_air_outlet_C, float(heater_state["airOutletTemp_C"]))
        max_wire_C = max(max_wire_C, float(heater_state["wireMaxTemp_C"]))
        max_hole_velocity_m_s = max(max_hole_velocity_m_s, float(branch_aeration["holeVelocity_m_s"]))
        max_flow_per_hole_lpm = max(max_flow_per_hole_lpm, float(branch_aeration["flowPerHole_Lpm"]))
        weighted_hwire += n_tubes_group * float(heater_state["hWireMean_W_m2K"])
        weighted_wire_mean += n_tubes_group * float(np.mean(heater_state["wireProfile_C"]))
        weighted_branch_q += n_tubes_group * float(branch_aeration["QtoBed_W"])
        tube_current_min_A = min(tube_current_min_A, tube_current_A)
        tube_current_max_A = max(tube_current_max_A, tube_current_A)
        tube_current_weighted_A += n_tubes_group * tube_current_A
        tube_power_min_W = min(tube_power_min_W, tube_power_W)
        tube_power_max_W = max(tube_power_max_W, tube_power_W)
        tube_power_weighted_W += n_tubes_group * tube_power_W

    mean_tube_current_A = tube_current_weighted_A / max(n_tubes, 1)
    mean_tube_power_W = tube_power_weighted_W / max(n_tubes, 1)
    mean_hwire = weighted_hwire / max(n_tubes, 1)
    mean_wire_temp_C = weighted_wire_mean / max(n_tubes, 1)
    mean_branch_q_W = weighted_branch_q / max(n_tubes, 1)
    equivalent_resistance_ohm = voltage_V / max(total_current_A, 1e-12)

    return {
        "nParallelTubes": int(n_tubes),
        "parallelStringCount": int(len(string_counts)),
        "stringTubeCounts": [int(count) for count in string_counts],
        "electricalTopologyLabel": topology_label,
        "branchFlow_Lpm": float(branch_flow_lpm),
        "branchResistance_Ohm": float(tube_resistance_ohm),
        "branchCurrent_A": float(mean_tube_current_A),
        "branchPower_W": float(mean_tube_power_W),
        "supplyEquivalentResistance_Ohm": float(equivalent_resistance_ohm),
        "tubeCurrentMin_A": float(tube_current_min_A),
        "tubeCurrentMean_A": float(mean_tube_current_A),
        "tubeCurrentMax_A": float(tube_current_max_A),
        "tubePowerMin_W": float(tube_power_min_W),
        "tubePowerMean_W": float(mean_tube_power_W),
        "tubePowerMax_W": float(tube_power_max_W),
        "hWire_W_m2K": float(mean_hwire),
        "totalPower_W": float(total_power_W),
        "totalCurrent_A": float(total_current_A),
        "branchQtoBed_W": float(mean_branch_q_W),
        "QtoBed_W": float(q_to_bed_W),
        "deltaP_Pa": float(max_branch_delta_p + distribution["total_Pa"]),
        "distributionDeltaP_Pa": float(distribution["total_Pa"]),
        "distributionHeader_Pa": float(distribution["header_Pa"]),
        "distributionHeaderToSplitter_Pa": float(distribution["header_to_splitter_Pa"]),
        "distributionSplitterBody_Pa": float(distribution["splitter_body_Pa"]),
        "distributionConnectorFriction_Pa": float(distribution["connector_friction_Pa"]),
        "distributionConnectorToBranch_Pa": float(distribution["connector_to_branch_Pa"]),
        "airOutlet_C": float(max_air_outlet_C),
        "wireMax_C": float(max_wire_C),
        "wireMean_C": float(mean_wire_temp_C),
        "flowPerHole_Lpm": float(max_flow_per_hole_lpm),
        "waterLoss_kg_day": float(water_loss_kg_day),
        "latentEvap_W": float(latent_evap_W),
        "holeVelocity_m_s": float(max_hole_velocity_m_s),
    }


IDELCHIK_CONICAL_ANGLE_DEG = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 60.0, 100.0, 140.0, 180.0], dtype=float)
IDELCHIK_CONICAL_EPSILON = np.array([0.025, 0.050, 0.075, 0.10, 0.15, 0.25, 0.60, 1.0], dtype=float)
IDELCHIK_CONICAL_C = np.array(
    [
        [1.00, 0.96, 0.93, 0.90, 0.86, 0.80, 0.69, 0.59, 0.50],
        [1.00, 0.93, 0.86, 0.80, 0.75, 0.67, 0.58, 0.53, 0.50],
        [1.00, 0.87, 0.75, 0.65, 0.58, 0.50, 0.48, 0.49, 0.50],
        [1.00, 0.80, 0.69, 0.55, 0.48, 0.41, 0.41, 0.44, 0.50],
        [1.00, 0.76, 0.58, 0.43, 0.33, 0.25, 0.27, 0.38, 0.50],
        [1.00, 0.68, 0.45, 0.30, 0.22, 0.17, 0.22, 0.34, 0.50],
        [1.00, 0.46, 0.27, 0.18, 0.14, 0.13, 0.21, 0.33, 0.50],
        [1.00, 0.32, 0.20, 0.14, 0.11, 0.10, 0.18, 0.30, 0.50],
    ],
    dtype=float,
)


def idelchik_conical_bellmouth_c(epsilon: float, cone_angle_deg: float) -> float:
    eps = max(min(float(epsilon), float(IDELCHIK_CONICAL_EPSILON[-1])), float(IDELCHIK_CONICAL_EPSILON[0]))
    angle = max(min(float(cone_angle_deg), float(IDELCHIK_CONICAL_ANGLE_DEG[-1])), float(IDELCHIK_CONICAL_ANGLE_DEG[0]))
    column_values = np.array(
        [np.interp(angle, IDELCHIK_CONICAL_ANGLE_DEG, row) for row in IDELCHIK_CONICAL_C],
        dtype=float,
    )
    return float(np.interp(eps, IDELCHIK_CONICAL_EPSILON, column_values))


def idelchik_conical_contraction_loss_Pa(
    mdot_kg_s: float,
    T_C: float,
    pressure_Pa: float,
    D_upstream_m: float,
    D_downstream_m: float,
    cone_angle_deg: float,
) -> float:
    if D_upstream_m <= 0.0 or D_downstream_m <= 0.0 or D_downstream_m >= D_upstream_m:
        return 0.0
    props = air_props(T_C, pressure_Pa)
    area_upstream_m2 = np.pi * D_upstream_m**2 / 4.0
    area_downstream_m2 = np.pi * D_downstream_m**2 / 4.0
    epsilon = area_downstream_m2 / max(area_upstream_m2, 1e-12)
    velocity_downstream_m_s = mdot_kg_s / max(props["rho_kg_m3"] * area_downstream_m2, 1e-12)
    # Idel'chik, Handbook of Hydraulic Resistance, Section III, Diagram 3-5:
    # conical converging bellmouth without end wall. The tabulated coefficient C
    # is used with Eq. (3-3): zeta = C * (1/epsilon - 1).
    coeff_c = idelchik_conical_bellmouth_c(epsilon, cone_angle_deg)
    loss_coeff = max(coeff_c * (1.0 / max(epsilon, 1e-12) - 1.0), 0.0)
    return float(loss_coeff * 0.5 * props["rho_kg_m3"] * velocity_downstream_m_s**2)


def sudden_expansion_loss_Pa(mdot_kg_s: float, T_C: float, pressure_Pa: float, D_upstream_m: float, D_downstream_m: float) -> float:
    if D_upstream_m <= 0.0 or D_downstream_m <= 0.0 or D_downstream_m <= D_upstream_m:
        return 0.0
    props = air_props(T_C, pressure_Pa)
    area_upstream_m2 = np.pi * D_upstream_m**2 / 4.0
    area_downstream_m2 = np.pi * D_downstream_m**2 / 4.0
    velocity_upstream_m_s = mdot_kg_s / max(props["rho_kg_m3"] * area_upstream_m2, 1e-12)
    epsilon = area_upstream_m2 / max(area_downstream_m2, 1e-12)
    # Idel'chik Section IV sudden expansion; identical to the standard
    # Borda-Carnot coefficient used in White/Munson when referenced to the
    # upstream velocity: zeta = (1 - A1/A2)^2.
    loss_coeff = max((1.0 - epsilon) ** 2, 0.0)
    return float(loss_coeff * 0.5 * props["rho_kg_m3"] * velocity_upstream_m_s**2)


def branch_distribution_losses(
    total_flow_m3_s: float,
    branch_flow_m3_s: float,
    T_C: float,
    pressure_Pa: float,
    n_tubes: int,
    aer_inputs: dict,
    downstream_branch_id_m: float,
) -> dict[str, float]:
    props = air_props(T_C, pressure_Pa)
    rho = props["rho_kg_m3"]
    mu = props["mu_Pa_s"]
    header_enabled = bool(aer_inputs.get("headerEnabled", True))
    header_minor_k = float(aer_inputs.get("headerMinorK", 1.0))
    header_id_m = max(1e-12, float(aer_inputs.get("headerID_mm", 50.0)) / 1000.0)
    splitter_outlet_count = max(1, int(round(float(aer_inputs.get("splitterOutletCount", 4)))))
    splitter_inlet_id_m = max(1e-12, float(aer_inputs.get("splitterInletID_mm", aer_inputs.get("branchConnectorID_mm", 19.0))) / 1000.0)
    branch_connector_id_m = max(1e-12, float(aer_inputs.get("branchConnectorID_mm", 19.0)) / 1000.0)
    branch_connector_length_m = max(0.0, float(aer_inputs.get("branchConnectorLength_m", 0.0)))
    splitter_body_k = max(0.0, float(aer_inputs.get("splitterBodyK", 1.0)))
    contraction_angle_deg = float(aer_inputs.get("contractionConeAngle_deg", 60.0))
    active_splitter_branches = max(1, min(splitter_outlet_count, int(n_tubes)))

    branch_mdot_kg_s = rho * branch_flow_m3_s
    splitter_inlet_flow_m3_s = active_splitter_branches * branch_flow_m3_s
    splitter_inlet_mdot_kg_s = rho * splitter_inlet_flow_m3_s

    connector_area_m2 = np.pi * branch_connector_id_m**2 / 4.0
    connector_velocity_m_s = branch_mdot_kg_s / max(rho * connector_area_m2, 1e-12)
    connector_re = rho * connector_velocity_m_s * branch_connector_id_m / max(mu, 1e-12)
    connector_f = churchill_friction_factor(connector_re)
    dp_connector_friction = (
        connector_f
        * branch_connector_length_m
        / max(branch_connector_id_m, 1e-12)
        * 0.5
        * rho
        * connector_velocity_m_s**2
    )

    if header_enabled:
        header_area_m2 = np.pi * header_id_m**2 / 4.0
        header_velocity_m_s = total_flow_m3_s / max(header_area_m2, 1e-12)
        dp_header = header_minor_k * 0.5 * rho * header_velocity_m_s**2
        dp_header_to_splitter = idelchik_conical_contraction_loss_Pa(
            splitter_inlet_mdot_kg_s,
            T_C,
            pressure_Pa,
            header_id_m,
            splitter_inlet_id_m,
            contraction_angle_deg,
        )
    else:
        dp_header = 0.0
        dp_header_to_splitter = 0.0

    dp_splitter_body = splitter_body_k * 0.5 * rho * connector_velocity_m_s**2
    dp_connector_to_branch = 0.0
    if downstream_branch_id_m < branch_connector_id_m:
        dp_connector_to_branch = idelchik_conical_contraction_loss_Pa(
            branch_mdot_kg_s,
            T_C,
            pressure_Pa,
            branch_connector_id_m,
            downstream_branch_id_m,
            contraction_angle_deg,
        )
    elif downstream_branch_id_m > branch_connector_id_m:
        dp_connector_to_branch = sudden_expansion_loss_Pa(
            branch_mdot_kg_s,
            T_C,
            pressure_Pa,
            branch_connector_id_m,
            downstream_branch_id_m,
        )

    return {
        "header_Pa": float(dp_header),
        "header_to_splitter_Pa": float(dp_header_to_splitter),
        "splitter_body_Pa": float(dp_splitter_body),
        "connector_friction_Pa": float(dp_connector_friction),
        "connector_to_branch_Pa": float(dp_connector_to_branch),
        "total_Pa": float(dp_header + dp_header_to_splitter + dp_splitter_body + dp_connector_friction + dp_connector_to_branch),
        "connector_velocity_m_s": float(connector_velocity_m_s),
        "splitter_active_branches": float(active_splitter_branches),
        "header_enabled": float(1.0 if header_enabled else 0.0),
    }


def build_bin_model_from_inputs(bin_inputs: dict, wall_inputs: dict, env_inputs: dict, design_bed_C: float) -> dict:
    L = float(bin_inputs.get("length_m", 1.22))
    W = float(bin_inputs.get("width_m", 0.61))
    H = float(bin_inputs.get("height_m", 0.61))
    fill_fraction = float(bin_inputs.get("fillFraction", 0.80))
    side_area_m2 = 2.0 * H * L + 2.0 * H * W
    bottom_area_m2 = L * W
    top_area_m2 = L * W
    total_area_m2 = side_area_m2 + bottom_area_m2 + top_area_m2
    vent_count = max(0, int(round(float(bin_inputs.get("ventHoleCount", 0)))))
    vent_diameter_m = max(0.0, float(bin_inputs.get("ventHoleDiameter_m", 0.0)))
    vent_hole_area_m2 = np.pi * vent_diameter_m**2 / 4.0
    vent_area_m2 = min(top_area_m2, vent_count * vent_hole_area_m2)
    covered_top_area_m2 = max(top_area_m2 - vent_area_m2, 0.0)

    h_ext = natural_convection_vertical_plate_h(
        design_bed_C,
        float(env_inputs.get("greenhouseAir_C", -15.0)),
        float(env_inputs.get("referencePlateLength_m", H)),
        float(env_inputs.get("pressure_Pa", 101325.0)),
    )
    u_wall = overall_wall_u(wall_inputs, h_ext)
    abottom_node_m2 = bottom_area_m2 + 0.5 * side_area_m2
    top_wall_area_m2 = 0.5 * side_area_m2

    if bool(bin_inputs.get("openTop", True)):
        top_char_length_m = top_area_m2 / max(2.0 * (L + W), 1e-12)
        h_top = natural_convection_horizontal_plate_h(
            design_bed_C,
            float(env_inputs.get("greenhouseAir_C", -15.0)),
            top_char_length_m,
            float(env_inputs.get("pressure_Pa", 101325.0)),
        )
        ua_top = u_wall * top_wall_area_m2 + h_top * top_area_m2
        ua_total = u_wall * (side_area_m2 + bottom_area_m2) + h_top * top_area_m2
        open_top_area_m2 = top_area_m2
        top_mode = "open top"
    elif vent_area_m2 > 0.0:
        vent_perimeter_m = max(np.pi * vent_diameter_m * vent_count, 1e-12)
        vent_char_length_m = vent_area_m2 / vent_perimeter_m
        h_top = natural_convection_horizontal_plate_h(
            design_bed_C,
            float(env_inputs.get("greenhouseAir_C", -15.0)),
            vent_char_length_m,
            float(env_inputs.get("pressure_Pa", 101325.0)),
        )
        ua_top = u_wall * (covered_top_area_m2 + top_wall_area_m2) + h_top * vent_area_m2
        ua_total = u_wall * (side_area_m2 + bottom_area_m2 + covered_top_area_m2) + h_top * vent_area_m2
        open_top_area_m2 = vent_area_m2
        top_mode = "covered top with vent hole"
    else:
        h_top = u_wall
        ua_top = u_wall * (top_area_m2 + top_wall_area_m2)
        ua_total = u_wall * total_area_m2
        open_top_area_m2 = 0.0
        top_mode = "covered top"

    volume_m3 = L * W * H * fill_fraction
    mass_kg = volume_m3 * float(bin_inputs.get("bulkDensity_kg_m3", 650.0))
    c_total = mass_kg * float(bin_inputs.get("cp_J_kgK", 3200.0))
    interface_area_m2 = L * W
    interface_thickness_m = max(H / 2.0, 1e-9)
    ua_internal = float(bin_inputs.get("k_W_mK", 0.45)) * interface_area_m2 / interface_thickness_m
    lc_m = volume_m3 / max(total_area_m2, 1e-12)
    bi_lumped = h_ext * lc_m / max(float(bin_inputs.get("k_W_mK", 0.45)), 1e-12)

    return {
        "totalArea_m2": total_area_m2,
        "topMode": top_mode,
        "topArea_m2": top_area_m2,
        "coveredTopArea_m2": covered_top_area_m2,
        "openTopArea_m2": open_top_area_m2,
        "ventHoleArea_m2": vent_area_m2,
        "externalH_W_m2K": h_ext,
        "topAmbientH_W_m2K": h_top,
        "Uwall_W_m2K": u_wall,
        "UAtotalAmbient_W_K": ua_total,
        "UAbottomAmbient_W_K": u_wall * abottom_node_m2,
        "UAtopAmbient_W_K": ua_top,
        "UAinternal_W_K": ua_internal,
        "Ctotal_J_K": c_total,
        "Cbottom_J_K": 0.5 * c_total,
        "Ctop_J_K": 0.5 * c_total,
        "mass_kg": mass_kg,
        "volume_m3": volume_m3,
        "characteristicLength_m": lc_m,
        "Bi_lumped": bi_lumped,
        "lumpedStrictlyValid": bi_lumped < 0.1,
    }


def lumped_loss_model_from_bin(bin_model: dict, ambient_C: float, design_bed_C: float) -> dict:
    ua = float(bin_model["UAtotalAmbient_W_K"])
    c_total = float(bin_model["Ctotal_J_K"])
    qreq = ua * (design_bed_C - ambient_C)
    tau_s = c_total / max(ua, 1e-12)
    return {
        "requiredHeat_W": qreq,
        "UA_W_K": ua,
        "tau_h": tau_s / 3600.0,
        "Bi_lumped": float(bin_model["Bi_lumped"]),
        "lumpedStrictlyValid": bool(bin_model["lumpedStrictlyValid"]),
        "externalH_W_m2K": float(bin_model["externalH_W_m2K"]),
    }


def two_node_steady_state_python(bin_model: dict, ambient_C: float, q_bed_W: float, bottom_fraction: float) -> tuple[float, float]:
    uab = float(bin_model["UAbottomAmbient_W_K"])
    uat = float(bin_model["UAtopAmbient_W_K"])
    u12 = float(bin_model["UAinternal_W_K"])
    qb = bottom_fraction * q_bed_W
    qt = (1.0 - bottom_fraction) * q_bed_W
    A = np.array([[uab + u12, -u12], [-u12, uat + u12]], dtype=float)
    b = np.array([qb + uab * ambient_C, qt + uat * ambient_C], dtype=float)
    x = np.linalg.solve(A, b)
    return float(x[0]), float(x[1])


def cooling_mode_config(config: dict) -> dict:
    summary_inputs = config.get("summary_inputs", {})
    model_overrides = normalized_model_overrides(config)
    flow_sweep = model_overrides.get("sweep", {})
    cooling = dict(config.get("cooling_mode", {}))
    cooling.setdefault("ambientAir_C", 41.0)
    cooling.setdefault("pressure_Pa", 101325.0)
    cooling.setdefault("designBed_C", float(summary_inputs.get("bedSetpoint_C", 20.0)))
    cooling.setdefault("method", "spot_cooler_with_assist_blower")
    cooling.setdefault("totalFlow_Lpm_min", float(flow_sweep.get("totalFlow_Lpm_min", 40.0)))
    cooling.setdefault("totalFlow_Lpm_max", float(flow_sweep.get("totalFlow_Lpm_max", 600.0)))
    cooling.setdefault("totalFlow_Lpm_count", int(flow_sweep.get("totalFlow_Lpm_count", 100)))
    cooling.setdefault(
        "inletRelativeHumidity",
        float(cooling.get("supplyRelativeHumidity", model_overrides.get("evaporation", {}).get("relativeHumidity", 0.8))),
    )
    spot = dict(cooling.get("spot_cooler", {}))
    spot.setdefault("targetSupplyAir_C", 18.0)
    spot.setdefault("maxCoolingCapacity_W", None)
    spot.setdefault("COP", None)
    cooling["spot_cooler"] = spot
    assist = dict(cooling.get("assist_blower", {}))
    assist.setdefault("enabled", True)
    assist.setdefault("model", "Whitewater WW-18 regenerative blower")
    assist.setdefault("ratedFlow_CFM", 8.0)
    assist.setdefault("ratedPressure_inH2O", 10.0)
    assist.setdefault("shutoffPressure_inH2O", 28.0)
    assist.setdefault("motorPower_W", 190.0)
    assist.setdefault(
        "reference",
        "AquaCave Whitewater WW-18 regenerative blower: 8 CFM at 10 in. water, max duty 28 in. water, 190 W",
    )
    cooling["assist_blower"] = assist
    opt = dict(cooling.get("optimization", {}))
    opt.setdefault("enabled", True)
    opt.setdefault("top_candidate_count", 12)
    opt.setdefault("n_top", 10)
    opt.setdefault("max_bottom_temp_C", float(model_overrides.get("bin", {}).get("maxSafe_C", 25.0)))
    opt.setdefault("max_top_temp_C", float(model_overrides.get("bin", {}).get("maxSafe_C", 25.0)))
    opt.setdefault("spread_limit_C", float(config.get("optimization", {}).get("spread_limit_C", 10.0)))
    opt.setdefault("priority_order", list(DEFAULT_COOLING_PRIORITY_ORDER))
    cooling["optimization"] = opt
    limits = dict(cooling.get("limits", {}))
    shared_limits = model_overrides.get("limits", {})
    limits.setdefault("maxPressureDrop_Pa", float(shared_limits.get("maxPressureDrop_Pa", 1500.0)))
    if "maxAssistBlowerDeliverablePressure_Pa" not in limits:
        if cooling["assist_blower"].get("enabled", False):
            limits["maxAssistBlowerDeliverablePressure_Pa"] = in_h2o_to_pa(
                float(cooling["assist_blower"].get("ratedPressure_inH2O", 10.0))
            )
        else:
            limits["maxAssistBlowerDeliverablePressure_Pa"] = np.inf
    if "maxWaterLoss_kg_day" not in limits:
        if "maxWaterSupply_kg_day" in limits:
            limits["maxWaterLoss_kg_day"] = float(limits["maxWaterSupply_kg_day"])
        else:
            limits["maxWaterLoss_kg_day"] = float(shared_limits.get("maxWaterLoss_kg_day", np.inf))
    limits.setdefault("maxHoleVelocity_m_s", float(cooling.get("maxHoleVelocity_m_s", np.inf)))
    cooling["limits"] = limits
    return cooling


def cost_analysis_config(config: dict) -> dict:
    cost = dict(config.get("cost_analysis", {}))
    cost.setdefault("enabled", True)
    cost.setdefault("electricity_price_usd_per_kWh", 0.2875)
    cost.setdefault("location", "Connecticut, USA")
    cost.setdefault("tariff_class", "residential")
    cost.setdefault(
        "source",
        "https://www.eia.gov/electricity/annual/html/epa_02_10.html ; Connecticut residential average retail electricity price for 2024 = 28.75 cents/kWh = $0.2875/kWh.",
    )
    return cost


def year_round_analysis_config(config: dict) -> dict:
    summary_inputs = config.get("summary_inputs", {})
    model_overrides = normalized_model_overrides(config)
    bin_inputs = effective_bin_inputs(summary_inputs, model_overrides)
    annual = dict(config.get("year_round_analysis", {}))
    heating_default_C = float(bin_inputs.get("setpoint_C", summary_inputs.get("bedSetpoint_C", 20.0)))
    cooling_default_C = float(cooling_mode_config(config).get("designBed_C", heating_default_C))
    annual.setdefault("enabled", True)
    annual.setdefault("location", "Connecticut, USA")
    annual.setdefault("months", ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
    annual["heatingSetpoint_C"] = (
        heating_default_C
        if annual.get("heatingSetpoint_C") in (None, "")
        else float(annual["heatingSetpoint_C"])
    )
    annual["coolingSetpoint_C"] = (
        cooling_default_C
        if annual.get("coolingSetpoint_C") in (None, "")
        else float(annual["coolingSetpoint_C"])
    )
    annual.setdefault("daily_profile_model", "periodic_gaussian_kernel_regression_with_event_overrides")
    climate = dict(annual.get("climate_data", {}))
    climate.setdefault("period", "2016-2025")
    climate.setdefault("stateCode", 6)
    climate.setdefault("stateAbbreviation", "CT")
    climate.setdefault("statewideMonthlyParameter", "tavg")
    climate.setdefault(
        "statewideMonthlySource",
        "https://www.ncei.noaa.gov/access/monitoring/climate-at-a-glance/statewide/time-series/6/tavg/1/{month}/2016-2025/data.json",
    )
    climate.setdefault(
        "statewideMonthlySourceDoc",
        "https://www.ncei.noaa.gov/monitoring-content/cag/text/time-series-api-doc.html",
    )
    climate.setdefault("eventStation", "USW00014740")
    climate.setdefault(
        "eventStationDescription",
        "NOAA daily-summaries station USW00014740 used as the Connecticut freeze/heat-wave timing proxy.",
    )
    climate.setdefault(
        "eventSource",
        "https://www.ncei.noaa.gov/access/services/data/v1?dataset=daily-summaries&stations=USW00014740&startDate=2016-01-01&endDate=2025-12-31&dataTypes=TMAX,TMIN,TAVG&units=standard&format=json",
    )
    climate.setdefault("daysPerMonth", [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    climate.setdefault("monthlyMean_C", [-1.34, 0.08, 3.89, 8.91, 14.78, 19.74, 23.32, 21.98, 18.16, 12.44, 5.54, 0.53])
    climate.setdefault("monthlyStd_C", [2.12, 1.31, 1.70, 1.39, 1.07, 0.78, 0.79, 1.22, 0.63, 1.09, 1.48, 1.66])
    climate.setdefault("baselineRegressionModel", "periodic_gaussian_kernel_regression")
    climate.setdefault("gaussianKernelBandwidth_days", 18.0)
    climate.setdefault(
        "baselineRegressionDescription",
        "Daily baseline ambient is reconstructed from the 12 monthly Connecticut mean-temperature anchors using periodic Gaussian kernel regression on day of year. Monthly standard deviations are regressed the same way for reporting and event-context purposes. Freeze and heat-wave overrides are then applied on top of that smooth baseline.",
    )
    freeze = dict(climate.get("freeze", {}))
    freeze.setdefault("definition", "Study freeze-day definition: NOAA daily-summaries TMIN <= 32 F.")
    freeze.setdefault("averageDaysPerMonth", [26.3, 24.1, 18.6, 6.3, 0.2, 0.0, 0.0, 0.0, 0.0, 2.3, 16.5, 25.2])
    freeze.setdefault("meanAnomalyBelowMonthlyMean_C", [0.79, 0.85, 1.84, 3.97, 6.45, 0.0, 0.0, 0.0, 0.0, 6.64, 2.51, 1.08])
    freeze.setdefault("meanDayOfMonth", [15.9, 13.9, 14.8, 10.1, 14.0, 0.0, 0.0, 0.0, 0.0, 24.2, 17.3, 16.1])
    climate["freeze"] = freeze
    heat_wave = dict(climate.get("heat_wave", {}))
    heat_wave.setdefault(
        "definition",
        "Study heat-wave definition: at least 3 consecutive days with NOAA daily-summaries TMAX >= 90 F.",
    )
    heat_wave.setdefault("averageSpellStartsPerMonth", [0.0, 0.0, 0.0, 0.0, 0.1, 0.7, 1.7, 0.9, 0.2, 0.0, 0.0, 0.0])
    heat_wave.setdefault("averageDaysPerMonth", [0.0, 0.0, 0.0, 0.0, 0.3, 2.8, 8.4, 3.9, 0.7, 0.0, 0.0, 0.0])
    heat_wave.setdefault("meanAnomalyAboveMonthlyMean_C", [0.0, 0.0, 0.0, 0.0, 10.96, 7.24, 3.83, 4.74, 8.43, 0.0, 0.0, 0.0])
    heat_wave.setdefault("meanStartDayOfMonth", [0.0, 0.0, 0.0, 0.0, 17.0, 18.7, 15.9, 11.3, 4.0, 0.0, 0.0, 0.0])
    climate["heat_wave"] = heat_wave
    annual["climate_data"] = climate
    annual["daysPerMonth"] = list(climate["daysPerMonth"])
    annual["monthlyAmbient_C"] = list(climate["monthlyMean_C"])
    annual.setdefault(
        "source",
        "Representative 365-day Connecticut climate year built from NOAA 2016-2025 statewide monthly mean temperatures and standard deviations, plus NOAA daily-station freeze and heat-wave statistics. The Python annual layer applies the current two-node bed model and re-solves the selected heating and cooling operating points at each day's ambient temperature before computing duty, electricity, water, and cost.",
    )
    return annual


def heating_reference_point(payload: dict) -> dict:
    ref = payload.get("recommended", {})
    if isinstance(ref, dict) and ref:
        return ref
    best = payload.get("best_available", {})
    return best if isinstance(best, dict) else {}


def cooling_reference_point(payload: dict) -> dict:
    ref = payload.get("recommended", {})
    if isinstance(ref, dict) and ref:
        return ref
    best = payload.get("best_available", {})
    return best if isinstance(best, dict) else {}


def largest_remainder_allocation(expected_values: list[float]) -> list[int]:
    expected = np.maximum(np.array([float(v) for v in expected_values], dtype=float), 0.0)
    base = np.floor(expected).astype(int)
    target_total = int(np.rint(np.sum(expected)))
    remainder = target_total - int(np.sum(base))
    if remainder > 0:
        order = np.argsort(-(expected - base))
        for idx in order[:remainder]:
            base[idx] += 1
    return base.tolist()


def gaussian_weighted_day_selection(
    n_days: int,
    count: int,
    center_day: float,
    spread_days: float,
    repulsion_days: float,
) -> list[int]:
    if n_days <= 0 or count <= 0:
        return []
    center = float(center_day) if float(center_day) > 0.0 else 0.5 * (n_days + 1)
    spread = max(float(spread_days), 1.0)
    repulsion = max(float(repulsion_days), 0.0)
    day_numbers = np.arange(1, n_days + 1, dtype=float)
    base_weights = np.exp(-0.5 * ((day_numbers - center) / spread) ** 2)
    selected: list[int] = []

    for _ in range(min(count, n_days)):
        scores = np.array(base_weights, copy=True)
        if selected and repulsion > 0.0:
            for chosen_day in selected:
                scores *= 1.0 - 0.65 * np.exp(-0.5 * ((day_numbers - chosen_day) / repulsion) ** 2)
        for chosen_day in selected:
            scores[int(chosen_day) - 1] = 0.0
        best_idx = int(np.argmax(scores))
        if float(scores[best_idx]) <= 0.0:
            break
        selected.append(best_idx + 1)

    return sorted(selected)


def centered_day_mask(n_days: int, count: int, center_day: float) -> np.ndarray:
    mask = np.zeros(max(n_days, 0), dtype=bool)
    spread_days = max(2.0, n_days / 6.0)
    repulsion_days = max(1.5, n_days / 14.0)
    for day in gaussian_weighted_day_selection(n_days, count, center_day, spread_days, repulsion_days):
        mask[day - 1] = True
    return mask


def heatwave_day_mask(n_days: int, spell_count: int, total_days: int, mean_start_day: float) -> tuple[np.ndarray, list[int], list[int]]:
    mask = np.zeros(max(n_days, 0), dtype=bool)
    if n_days <= 0 or spell_count <= 0:
        return mask, [], []

    max_spells = max((n_days + 1) // 4, 0)
    spell_count = min(spell_count, max_spells)
    if spell_count <= 0:
        return mask, [], []

    min_required_days = 3 * spell_count
    max_assignable_days = max(n_days - (spell_count - 1), min_required_days)
    total_days = int(max(total_days, min_required_days))
    total_days = min(total_days, max_assignable_days)

    lengths = [3] * spell_count
    extra_days = total_days - 3 * spell_count
    idx = 0
    while extra_days > 0:
        lengths[idx % spell_count] += 1
        extra_days -= 1
        idx += 1

    center = float(mean_start_day) if float(mean_start_day) > 0.0 else 0.5 * (n_days + 1)
    spacing = max(4.0, float(np.mean(lengths)) + 1.0)
    raw_starts = gaussian_weighted_day_selection(
        n_days,
        spell_count,
        center,
        max(2.0, n_days / 7.0),
        spacing,
    )
    if len(raw_starts) < spell_count:
        fallback = sorted(set(raw_starts + [int(round(center))]))
        while len(fallback) < spell_count:
            candidate = min(max(len(fallback) + 1, 1), n_days)
            if candidate not in fallback:
                fallback.append(candidate)
            else:
                break
        raw_starts = fallback[:spell_count]

    starts: list[int] = []
    for spell_idx, length in enumerate(lengths):
        earliest = 1 if not starts else starts[-1] + lengths[spell_idx - 1] + 1
        remaining_lengths = sum(lengths[spell_idx:])
        remaining_gaps = spell_count - spell_idx - 1
        latest = n_days - (remaining_lengths + remaining_gaps) + 1
        latest = max(earliest, latest)
        proposed = raw_starts[spell_idx] if spell_idx < len(raw_starts) else earliest
        starts.append(int(np.clip(proposed, earliest, latest)))

    for start, length in zip(starts, lengths):
        end = min(n_days, start + length - 1)
        mask[start - 1 : end] = True

    return mask, starts, lengths


def circular_gaussian_kernel_regression(
    x_eval: np.ndarray,
    x_anchor: np.ndarray,
    y_anchor: np.ndarray,
    bandwidth_days: float,
    period_days: float,
) -> np.ndarray:
    bandwidth = max(float(bandwidth_days), 1e-6)
    result = np.zeros_like(x_eval, dtype=float)
    for idx, x_val in enumerate(x_eval):
        deltas = np.abs(x_anchor - x_val)
        deltas = np.minimum(deltas, period_days - deltas)
        weights = np.exp(-0.5 * (deltas / bandwidth) ** 2)
        result[idx] = float(np.dot(weights, y_anchor) / max(np.sum(weights), 1e-12))
    return result


def build_representative_climate_year(annual_cfg: dict) -> dict:
    climate = annual_cfg["climate_data"]
    months = list(annual_cfg.get("months", []))
    days_per_month = [int(v) for v in climate.get("daysPerMonth", [])]
    monthly_mean = [float(v) for v in climate.get("monthlyMean_C", [])]
    monthly_std = [float(v) for v in climate.get("monthlyStd_C", [])]
    freeze_cfg = dict(climate.get("freeze", {}))
    heat_cfg = dict(climate.get("heat_wave", {}))
    freeze_expected = [float(v) for v in freeze_cfg.get("averageDaysPerMonth", [])]
    freeze_anomaly = [float(v) for v in freeze_cfg.get("meanAnomalyBelowMonthlyMean_C", [])]
    freeze_mean_day = [float(v) for v in freeze_cfg.get("meanDayOfMonth", [])]
    heat_expected_starts = [float(v) for v in heat_cfg.get("averageSpellStartsPerMonth", [])]
    heat_expected_days = [float(v) for v in heat_cfg.get("averageDaysPerMonth", [])]
    heat_anomaly = [float(v) for v in heat_cfg.get("meanAnomalyAboveMonthlyMean_C", [])]
    heat_mean_start = [float(v) for v in heat_cfg.get("meanStartDayOfMonth", [])]

    n_months = min(
        len(months),
        len(days_per_month),
        len(monthly_mean),
        len(monthly_std),
        len(freeze_expected),
        len(freeze_anomaly),
        len(freeze_mean_day),
        len(heat_expected_starts),
        len(heat_expected_days),
        len(heat_anomaly),
        len(heat_mean_start),
    )
    months = months[:n_months]
    days_per_month = days_per_month[:n_months]
    monthly_mean = monthly_mean[:n_months]
    monthly_std = monthly_std[:n_months]
    freeze_expected = freeze_expected[:n_months]
    freeze_anomaly = freeze_anomaly[:n_months]
    freeze_mean_day = freeze_mean_day[:n_months]
    heat_expected_starts = heat_expected_starts[:n_months]
    heat_expected_days = heat_expected_days[:n_months]
    heat_anomaly = heat_anomaly[:n_months]
    heat_mean_start = heat_mean_start[:n_months]

    freeze_alloc = largest_remainder_allocation(freeze_expected)
    heat_start_alloc = largest_remainder_allocation(heat_expected_starts)
    heat_day_alloc = largest_remainder_allocation(heat_expected_days)
    total_days = int(sum(days_per_month))
    month_midpoints = []
    cursor = 0.0
    for n_days in days_per_month:
        month_midpoints.append(cursor + 0.5 * (n_days + 1))
        cursor += float(n_days)
    month_midpoints_arr = np.array(month_midpoints, dtype=float)
    mean_arr = np.array(monthly_mean, dtype=float)
    std_arr = np.array(monthly_std, dtype=float)
    bandwidth_days = float(climate.get("gaussianKernelBandwidth_days", 18.0))
    day_positions = np.arange(1, total_days + 1, dtype=float)
    daily_mean_regressed = circular_gaussian_kernel_regression(
        day_positions,
        month_midpoints_arr,
        mean_arr,
        bandwidth_days,
        float(total_days),
    )
    daily_std_regressed = circular_gaussian_kernel_regression(
        day_positions,
        month_midpoints_arr,
        std_arr,
        bandwidth_days,
        float(total_days),
    )

    daily_rows: list[dict] = []
    monthly_rows: list[dict] = []
    day_of_year = 0

    for month_idx, month in enumerate(months):
        n_days = int(days_per_month[month_idx])
        mean_c = float(monthly_mean[month_idx])
        std_c = float(monthly_std[month_idx])
        day_slice = slice(day_of_year, day_of_year + n_days)
        base_profile = np.array(daily_mean_regressed[day_slice], copy=True)
        sigma_profile = np.array(daily_std_regressed[day_slice], copy=True)

        freeze_mask = centered_day_mask(n_days, int(freeze_alloc[month_idx]), float(freeze_mean_day[month_idx]))
        heat_mask, heat_starts, heat_lengths = heatwave_day_mask(
            n_days,
            int(heat_start_alloc[month_idx]),
            int(heat_day_alloc[month_idx]),
            float(heat_mean_start[month_idx]),
        )

        ambient_profile = np.array(base_profile, copy=True)
        if float(freeze_anomaly[month_idx]) > 0.0 and np.any(freeze_mask):
            ambient_profile[freeze_mask] = np.minimum(
                ambient_profile[freeze_mask],
                mean_c - float(freeze_anomaly[month_idx]),
            )
        if float(heat_anomaly[month_idx]) > 0.0 and np.any(heat_mask):
            ambient_profile[heat_mask] = np.maximum(
                ambient_profile[heat_mask],
                mean_c + float(heat_anomaly[month_idx]),
            )

        monthly_rows.append(
            {
                "month": month,
                "monthIndex": month_idx,
                "days": n_days,
                "monthlyMean_C": mean_c,
                "monthlyStd_C": std_c,
                "gaussianBandwidth_days": bandwidth_days,
                "freezeDaysExpected": float(freeze_expected[month_idx]),
                "freezeDaysAllocated": int(np.count_nonzero(freeze_mask)),
                "freezeMeanDayOfMonth": float(freeze_mean_day[month_idx]),
                "freezeAnomalyBelowMean_C": float(freeze_anomaly[month_idx]),
                "heatWaveSpellStartsExpected": float(heat_expected_starts[month_idx]),
                "heatWaveSpellStartsAllocated": len(heat_starts),
                "heatWaveDaysExpected": float(heat_expected_days[month_idx]),
                "heatWaveDaysAllocated": int(np.count_nonzero(heat_mask)),
                "heatWaveMeanStartDayOfMonth": float(heat_mean_start[month_idx]),
                "heatWaveAnomalyAboveMean_C": float(heat_anomaly[month_idx]),
                "heatWaveStartDays": [int(v) for v in heat_starts],
                "heatWaveSpellLengths_days": [int(v) for v in heat_lengths],
                "dayOfYearStart": day_of_year + 1,
                "dayOfYearEnd": day_of_year + n_days,
                "dayOfYearMid": day_of_year + 0.5 * (n_days + 1),
            }
        )

        for local_idx in range(n_days):
            day_of_year += 1
            daily_rows.append(
                {
                    "dayOfYear": int(day_of_year),
                    "monthIndex": month_idx,
                    "month": month,
                    "dayOfMonth": int(local_idx + 1),
                    "ambientBase_C": float(base_profile[local_idx]),
                    "ambientSigma_C": float(sigma_profile[local_idx]),
                    "ambient_C": float(ambient_profile[local_idx]),
                    "freezeDay": bool(freeze_mask[local_idx]),
                    "heatWaveDay": bool(heat_mask[local_idx]),
                    "monthlyMean_C": mean_c,
                    "monthlyStd_C": std_c,
                }
            )

    return {
        "daily": daily_rows,
        "monthly": monthly_rows,
        "annualFreezeDays": int(sum(item["freezeDaysAllocated"] for item in monthly_rows)),
        "annualHeatWaveDays": int(sum(item["heatWaveDaysAllocated"] for item in monthly_rows)),
        "annualHeatWaveSpellStarts": int(sum(item["heatWaveSpellStartsAllocated"] for item in monthly_rows)),
    }


def simulate_year_round_day_task(
    task: tuple[int, dict, dict, dict, dict, float, float, float]
) -> tuple[int, dict]:
    idx, climate_day, config, heating_ref, cooling_ref, heating_setpoint_C, cooling_setpoint_C, electricity_price = task
    summary_inputs = config.get("summary_inputs", {})
    model_overrides = normalized_model_overrides(config)
    env_inputs_base = effective_environment_inputs(summary_inputs, model_overrides)
    bin_inputs = effective_bin_inputs(summary_inputs, model_overrides)
    wall_inputs = effective_bin_wall_inputs(summary_inputs, model_overrides)
    bottom_fraction = float(bin_inputs.get("heatToBottomFraction", 0.60))

    ambient_C = float(climate_day["ambient_C"])
    env_inputs = dict(env_inputs_base)
    env_inputs["greenhouseAir_C"] = ambient_C
    bin_model_heat = build_bin_model_from_inputs(bin_inputs, wall_inputs, env_inputs, heating_setpoint_C)
    bin_model_cool = build_bin_model_from_inputs(bin_inputs, wall_inputs, env_inputs, cooling_setpoint_C)
    required_heat_W = max(float(lumped_loss_model_from_bin(bin_model_heat, ambient_C, heating_setpoint_C)["requiredHeat_W"]), 0.0)
    required_cooling_W = max(-float(lumped_loss_model_from_bin(bin_model_cool, ambient_C, cooling_setpoint_C)["requiredHeat_W"]), 0.0)

    mode = "idle"
    nominal_capacity_W = 0.0
    nominal_power_W = 0.0
    water_rate_kg_day = 0.0
    required_load_W = 0.0
    nominal_state: dict[str, float] = {}
    if required_heat_W > 0.0:
        mode = "heating"
        nominal_state = simulate_heating_reference_point(config, heating_ref, ambient_C, heating_setpoint_C)
        nominal_capacity_W = max(float(nominal_state.get("QtoBed_W", 0.0)), 0.0)
        nominal_power_W = max(float(nominal_state.get("totalPower_W", 0.0)), 0.0)
        water_rate_kg_day = max(float(nominal_state.get("waterLoss_kg_day", 0.0)), 0.0)
        required_load_W = required_heat_W
    elif required_cooling_W > 0.0:
        mode = "cooling"
        nominal_state = simulate_summer_cooling_point(
            config,
            float(cooling_ref.get("totalFlow_Lpm", 0.0)),
            bin_model_cool,
            required_cooling_W,
            ambient_air_override_C=ambient_C,
        )
        nominal_capacity_W = max(-float(nominal_state.get("QtoBed_W", 0.0)), 0.0)
        nominal_power_W = max(
            float(nominal_state.get("spotCoolerPower_W", 0.0)) + float(nominal_state.get("assistBlowerPower_W", 0.0)),
            0.0,
        )
        water_rate_kg_day = max(float(nominal_state.get("waterLoss_kg_day", 0.0)), 0.0)
        required_load_W = required_cooling_W

    duty_required = required_load_W / max(nominal_capacity_W, 1e-12) if nominal_capacity_W > 0.0 else (0.0 if required_load_W <= 0.0 else np.inf)
    duty_applied = min(max(duty_required, 0.0), 1.0) if np.isfinite(duty_required) else 1.0
    unmet_load_W = max(required_load_W - nominal_capacity_W, 0.0) if required_load_W > 0.0 else 0.0

    if mode == "heating":
        q_bed_applied_W = nominal_capacity_W * duty_applied
        tb_C, tt_C = two_node_steady_state_python(bin_model_heat, ambient_C, q_bed_applied_W, bottom_fraction)
        heating_kWh = nominal_power_W * duty_applied * 24.0 / 1000.0
        cooling_kWh = 0.0
    elif mode == "cooling":
        q_bed_applied_W = -nominal_capacity_W * duty_applied
        tb_C, tt_C = two_node_steady_state_python(bin_model_cool, ambient_C, q_bed_applied_W, bottom_fraction)
        heating_kWh = 0.0
        cooling_kWh = nominal_power_W * duty_applied * 24.0 / 1000.0
    else:
        q_bed_applied_W = 0.0
        tb_C, tt_C = two_node_steady_state_python(bin_model_heat, ambient_C, 0.0, bottom_fraction)
        heating_kWh = 0.0
        cooling_kWh = 0.0

    total_kWh = heating_kWh + cooling_kWh
    cost_usd = total_kWh * electricity_price
    water_kg = water_rate_kg_day * duty_applied
    unmet_kWh_equivalent = unmet_load_W * 24.0 / 1000.0

    day_record = {
        "dayOfYear": int(climate_day["dayOfYear"]),
        "monthIndex": int(climate_day["monthIndex"]),
        "month": str(climate_day["month"]),
        "dayOfMonth": int(climate_day["dayOfMonth"]),
        "ambientBase_C": float(climate_day["ambientBase_C"]),
        "ambientSigma_C": float(climate_day["ambientSigma_C"]),
        "ambient_C": ambient_C,
        "freezeDay": bool(climate_day["freezeDay"]),
        "heatWaveDay": bool(climate_day["heatWaveDay"]),
        "monthlyMean_C": float(climate_day["monthlyMean_C"]),
        "monthlyStd_C": float(climate_day["monthlyStd_C"]),
        "mode": mode,
        "requiredHeat_W": float(required_heat_W),
        "requiredCooling_W": float(required_cooling_W),
        "requiredLoad_W": float(required_load_W),
        "nominalCapacity_W": float(nominal_capacity_W),
        "nominalPower_W": float(nominal_power_W),
        "dutyRequired": float(duty_required),
        "dutyApplied": float(duty_applied),
        "unmetLoad_W": float(unmet_load_W),
        "unmet_kWh_equivalent": float(unmet_kWh_equivalent),
        "qBedApplied_W": float(q_bed_applied_W),
        "bottom_C": float(tb_C),
        "top_C": float(tt_C),
        "heating_kWh": float(heating_kWh),
        "cooling_kWh": float(cooling_kWh),
        "total_kWh": float(total_kWh),
        "cost_usd": float(cost_usd),
        "water_kg": float(water_kg),
        "nominalState": nominal_state,
    }
    return idx, day_record


def build_year_round_payload(heating_payload: dict, cooling_payload: dict, config: dict) -> dict:
    annual_cfg = year_round_analysis_config(config)
    cost_cfg = cost_analysis_config(config)
    parallel_cfg = python_parallel_config(config)
    climate_payload = build_representative_climate_year(annual_cfg)

    heating_ref = heating_reference_point(heating_payload)
    cooling_ref = cooling_reference_point(cooling_payload)
    heating_setpoint_C = float(annual_cfg["heatingSetpoint_C"])
    cooling_setpoint_C = float(annual_cfg["coolingSetpoint_C"])
    electricity_price = float(cost_cfg["electricity_price_usd_per_kWh"])

    daily: list[dict] = []
    monthly: list[dict] = []
    annual_heating_kWh = 0.0
    annual_cooling_kWh = 0.0
    annual_cost_usd = 0.0
    annual_water_kg = 0.0
    annual_unmet_kWh = 0.0
    day_tasks = [
        (
            idx,
            climate_day,
            config,
            heating_ref,
            cooling_ref,
            heating_setpoint_C,
            cooling_setpoint_C,
            electricity_price,
        )
        for idx, climate_day in enumerate(climate_payload["daily"])
    ]
    use_parallel = bool(parallel_cfg.get("enabled", False) and parallel_cfg.get("yearRound", True) and len(day_tasks) > 1)
    if use_parallel:
        ctx = mp.get_context("spawn")
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=python_process_pool_workers(parallel_cfg),
            mp_context=ctx,
        ) as executor:
            day_results = list(executor.map(simulate_year_round_day_task, day_tasks))
    else:
        day_results = [simulate_year_round_day_task(task) for task in day_tasks]

    day_results.sort(key=lambda item: item[0])
    for _, day_record in day_results:
        annual_heating_kWh += float(day_record["heating_kWh"])
        annual_cooling_kWh += float(day_record["cooling_kWh"])
        annual_cost_usd += float(day_record["cost_usd"])
        annual_water_kg += float(day_record["water_kg"])
        annual_unmet_kWh += float(day_record["unmet_kWh_equivalent"])
        day_record["appliedQtoBed_W"] = float(day_record.pop("qBedApplied_W"))
        day_record["nominalAirOutlet_C"] = float(day_record["nominalState"].get("airOutlet_C", np.nan)) if day_record["nominalState"] else float("nan")
        day_record["nominalWireMax_C"] = float(day_record["nominalState"].get("wireMax_C", np.nan)) if day_record["nominalState"] else float("nan")
        day_record["nominalWaterLoss_kg_day"] = float(day_record["nominalState"].get("waterLoss_kg_day", 0.0)) if day_record["nominalState"] else 0.0
        daily.append(day_record)

    cumulative_cost = 0.0
    cumulative_energy = 0.0
    for item in daily:
        cumulative_cost += float(item["cost_usd"])
        cumulative_energy += float(item["total_kWh"])
        item["cumulativeCost_usd"] = float(cumulative_cost)
        item["cumulativeEnergy_kWh"] = float(cumulative_energy)

    for climate_month in climate_payload["monthly"]:
        month_days = [item for item in daily if int(item["monthIndex"]) == int(climate_month["monthIndex"])]
        if not month_days:
            continue
        ambient_vals = np.array([float(item["ambient_C"]) for item in month_days], dtype=float)
        bottom_vals = np.array([float(item["bottom_C"]) for item in month_days], dtype=float)
        top_vals = np.array([float(item["top_C"]) for item in month_days], dtype=float)
        duty_vals = np.array([float(item["dutyApplied"]) for item in month_days], dtype=float)
        monthly.append(
            {
                "month": climate_month["month"],
                "days": int(climate_month["days"]),
                "ambientMean_C": float(np.mean(ambient_vals)),
                "ambientMin_C": float(np.min(ambient_vals)),
                "ambientMax_C": float(np.max(ambient_vals)),
                "monthlyMean_C": float(climate_month["monthlyMean_C"]),
                "monthlyStd_C": float(climate_month["monthlyStd_C"]),
                "freezeDays": int(sum(1 for item in month_days if item["freezeDay"])),
                "heatWaveDays": int(sum(1 for item in month_days if item["heatWaveDay"])),
                "heatWaveSpellStarts": int(climate_month["heatWaveSpellStartsAllocated"]),
                "bottomMean_C": float(np.mean(bottom_vals)),
                "topMean_C": float(np.mean(top_vals)),
                "bottomMin_C": float(np.min(bottom_vals)),
                "topMax_C": float(np.max(top_vals)),
                "dutyMean": float(np.mean(duty_vals)),
                "dutyPeak": float(np.max(duty_vals)),
                "heating_kWh": float(sum(float(item["heating_kWh"]) for item in month_days)),
                "cooling_kWh": float(sum(float(item["cooling_kWh"]) for item in month_days)),
                "total_kWh": float(sum(float(item["total_kWh"]) for item in month_days)),
                "cost_usd": float(sum(float(item["cost_usd"]) for item in month_days)),
                "water_kg": float(sum(float(item["water_kg"]) for item in month_days)),
                "unmetPeak_W": float(max(float(item["unmetLoad_W"]) for item in month_days)),
            }
        )

    return {
        "daily": daily,
        "months": monthly,
        "climate": climate_payload,
        "annualHeating_kWh": float(annual_heating_kWh),
        "annualCooling_kWh": float(annual_cooling_kWh),
        "annualTotal_kWh": float(annual_heating_kWh + annual_cooling_kWh),
        "annualCost_usd": float(annual_cost_usd),
        "annualWater_kg": float(annual_water_kg),
        "annualUnmet_kWh_equivalent": float(annual_unmet_kWh),
        "electricityPrice_usd_per_kWh": electricity_price,
        "costConfig": cost_cfg,
        "yearRoundConfig": annual_cfg,
        "heatingReference": heating_ref,
        "coolingReference": cooling_ref,
    }


def plot_year_round_energy_cost(payload: dict, output_dir: Path) -> None:
    daily = payload["daily"]
    climate_months = payload["climate"]["monthly"]
    x = np.array([int(item["dayOfYear"]) for item in daily], dtype=float)
    ambient = np.array([float(item["ambient_C"]) for item in daily], dtype=float)
    ambient_base = np.array([float(item["ambientBase_C"]) for item in daily], dtype=float)
    bottom = np.array([float(item["bottom_C"]) for item in daily], dtype=float)
    top = np.array([float(item["top_C"]) for item in daily], dtype=float)
    daily_cost = np.array([float(item["cost_usd"]) for item in daily], dtype=float)
    cumulative_cost = np.array([float(item["cumulativeCost_usd"]) for item in daily], dtype=float)
    freeze_mask = np.array([bool(item["freezeDay"]) for item in daily], dtype=bool)
    heat_mask = np.array([bool(item["heatWaveDay"]) for item in daily], dtype=bool)
    heating_mask = np.array([str(item["mode"]) == "heating" for item in daily], dtype=bool)
    cooling_mask = np.array([str(item["mode"]) == "cooling" for item in daily], dtype=bool)

    fig, axes = plt.subplots(2, 1, figsize=(14, 10), constrained_layout=True)

    ax = axes[0]
    ax.plot(x, ambient_base, color="#9da4aa", linewidth=1.1, linestyle="--", label="Gaussian-regressed ambient baseline")
    ax.plot(x, ambient, color="tab:gray", linewidth=1.6, label="Ambient with event overrides")
    ax.plot(x, bottom, color="tab:blue", linewidth=1.9, label="Bottom node")
    ax.plot(x, top, color="tab:red", linewidth=1.9, label="Top node")
    ax.axhspan(15.0, 25.0, color="#dfead7", alpha=0.40, label="Worm-safe band")
    if np.any(freeze_mask):
        ax.scatter(x[freeze_mask], ambient[freeze_mask], marker="v", s=12, color="#1f77b4", label="Freeze days")
    if np.any(heat_mask):
        ax.scatter(x[heat_mask], ambient[heat_mask], marker="^", s=14, color="#d62728", label="Heat-wave days")
    for month_info in climate_months[:-1]:
        ax.axvline(float(month_info["dayOfYearEnd"]) + 0.5, color="#d0d0d0", linewidth=0.8, alpha=0.7)
    ax.set_xticks([float(item["dayOfYearMid"]) for item in climate_months], [item["month"] for item in climate_months])
    ax.set_ylabel("Temperature (C)")
    ax.set_title("Representative Connecticut Climate Year and Daily Bed Temperatures")
    ax.grid(True, alpha=0.22)
    ax.legend(ncol=3, fontsize=8.5, loc="upper right")

    ax = axes[1]
    bar_colors = np.where(heating_mask, "#c66b4e", np.where(cooling_mask, "#4f8cc9", "#b7b7b7"))
    ax.bar(x, daily_cost, width=1.0, color=bar_colors, alpha=0.85)
    ax.set_ylabel("Daily cost (USD/day)")
    ax.set_xticks([float(item["dayOfYearMid"]) for item in climate_months], [item["month"] for item in climate_months])
    ax.grid(True, axis="y", alpha=0.22)
    for month_info in climate_months[:-1]:
        ax.axvline(float(month_info["dayOfYearEnd"]) + 0.5, color="#d0d0d0", linewidth=0.8, alpha=0.7)
    ax2 = ax.twinx()
    ax2.plot(x, cumulative_cost, color="tab:green", linewidth=2.0, label="Cumulative cost")
    ax2.set_ylabel("Cumulative cost (USD)")
    handles1 = [
        Patch(facecolor="#c66b4e", edgecolor="#c66b4e", alpha=0.85, label="Daily heating cost"),
        Patch(facecolor="#4f8cc9", edgecolor="#4f8cc9", alpha=0.85, label="Daily cooling cost"),
        Patch(facecolor="#b7b7b7", edgecolor="#b7b7b7", alpha=0.85, label="Daily idle cost"),
    ]
    labels1 = [handle.get_label() for handle in handles1]
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(handles1 + handles2, labels1 + labels2, fontsize=9, loc="upper left")
    ax.set_title("Daily and Cumulative Operating Cost vs Time")

    fig.suptitle(
        "Year-Round Temperature and Cost vs Time\n"
        f"JSON inputs used: T_heat,set = {float(payload['yearRoundConfig']['heatingSetpoint_C']):.1f} C, "
        f"T_cool,set = {float(payload['yearRoundConfig']['coolingSetpoint_C']):.1f} C, "
        f"electricity = {float(payload['electricityPrice_usd_per_kWh']):.4f} USD/kWh",
        fontsize=15,
        fontweight="bold",
    )
    fig.savefig(output_dir / "year_round_energy_cost_analysis.png", dpi=220)
    plt.close(fig)


def write_year_round_summary(payload: dict, output_dir: Path) -> None:
    cost_cfg = payload["costConfig"]
    annual_cfg = payload["yearRoundConfig"]
    climate_cfg = annual_cfg["climate_data"]
    freeze_cfg = climate_cfg["freeze"]
    heat_cfg = climate_cfg["heat_wave"]
    lines = [
        "============================================================",
        "YEAR-ROUND ENERGY, COST, AND WATER ANALYSIS",
        "============================================================",
        f"Location                             : {annual_cfg['location']}",
        f"Electricity price basis             : {float(payload['electricityPrice_usd_per_kWh']):.4f} USD/kWh ({cost_cfg['tariff_class']}, {cost_cfg['location']})",
        f"Electricity price source            : {cost_cfg['source']}",
        f"Climate baseline                    : {annual_cfg['source']}",
        f"NOAA statewide monthly source       : {climate_cfg['statewideMonthlySource']}",
        f"NOAA monthly API documentation      : {climate_cfg['statewideMonthlySourceDoc']}",
        f"NOAA daily event source             : {climate_cfg['eventSource']}",
        f"Representative period               : {climate_cfg['period']}",
        f"Representative profile model        : {annual_cfg['daily_profile_model']}",
        f"Resolved heating setpoint           : {float(annual_cfg['heatingSetpoint_C']):.1f} C",
        f"Resolved cooling setpoint           : {float(annual_cfg['coolingSetpoint_C']):.1f} C",
        f"Baseline regression model           : {climate_cfg['baselineRegressionModel']}",
        f"Gaussian kernel bandwidth           : {float(climate_cfg['gaussianKernelBandwidth_days']):.1f} days",
        f"Baseline regression description     : {climate_cfg['baselineRegressionDescription']}",
        f"Freeze-day definition               : {freeze_cfg['definition']}",
        f"Heat-wave definition                : {heat_cfg['definition']}",
        f"Allocated annual freeze days        : {int(payload['climate']['annualFreezeDays'])}",
        f"Allocated annual heat-wave starts   : {int(payload['climate']['annualHeatWaveSpellStarts'])}",
        f"Allocated annual heat-wave days     : {int(payload['climate']['annualHeatWaveDays'])}",
        "Daily power evaluation              : Heating and cooling electricity are computed from ambient-specific daily reevaluation of the selected operating points, not from fixed -15 C or 37.5 C reference power values.",
        "",
        "Reference operating points:",
        f"  Heating reference (design point)  : {float(payload['heatingReference'].get('totalPower_W', np.nan)):.1f} W electric, {float(payload['heatingReference'].get('QtoBed_W', np.nan)):.1f} W to bed",
        f"  Cooling reference (design point)  : {float(payload['coolingReference'].get('spotCoolerPower_W', 0.0) + payload['coolingReference'].get('assistBlowerPower_W', 0.0)):.1f} W electric, {max(-float(payload['coolingReference'].get('QtoBed_W', 0.0)), 0.0):.1f} W cooling to bed",
        "",
        "Annual totals:",
        f"  Heating electricity               : {float(payload['annualHeating_kWh']):.1f} kWh",
        f"  Cooling electricity               : {float(payload['annualCooling_kWh']):.1f} kWh",
        f"  Total electricity                 : {float(payload['annualTotal_kWh']):.1f} kWh",
        f"  Annual operating cost             : {float(payload['annualCost_usd']):.2f} USD",
        f"  Annual water replacement          : {float(payload['annualWater_kg']):.1f} kg",
        f"  Annual unmet-load equivalent      : {float(payload['annualUnmet_kWh_equivalent']):.1f} kWh-equivalent",
        "",
        "Monthly climate/event table:",
        "  Month  Tmean(C)  Tstd(C)  Tamb,min  Tamb,max  Freeze(d)  HW starts  HW days",
    ]

    for item in payload["months"]:
        lines.append(
            f"  {item['month']:>3}   {float(item['monthlyMean_C']):7.2f}  {float(item['monthlyStd_C']):6.2f}  "
            f"{float(item['ambientMin_C']):8.2f}  {float(item['ambientMax_C']):8.2f}  "
            f"{int(item['freezeDays']):9d}  {int(item['heatWaveSpellStarts']):9d}  {int(item['heatWaveDays']):7d}"
        )

    lines.extend(
        [
            "",
            "Monthly operating table:",
            "  Month  Tbottom,mean  Ttop,mean  Tbottom,min  Ttop,max  Duty,mean  Duty,peak  Eheat(kWh)  Ecool(kWh)  Cost(USD)  Water(kg)  Unmet,peak(W)",
        ]
    )

    for item in payload["months"]:
        lines.append(
            f"  {item['month']:>3}   {float(item['bottomMean_C']):11.2f}  {float(item['topMean_C']):9.2f}  "
            f"{float(item['bottomMin_C']):11.2f}  {float(item['topMax_C']):8.2f}  "
            f"{float(item['dutyMean']):9.3f}  {float(item['dutyPeak']):9.3f}  "
            f"{float(item['heating_kWh']):10.1f}  {float(item['cooling_kWh']):10.1f}  "
            f"{float(item['cost_usd']):9.2f}  {float(item['water_kg']):9.2f}  {float(item['unmetPeak_W']):12.1f}"
        )

    lines.append("============================================================")
    (output_dir / "study_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def interpolate_profile_value(distances_m: np.ndarray, values: np.ndarray, distance_m: float) -> float:
    if values.size == 0:
        return float("nan")
    if distance_m <= float(distances_m[0]):
        return float(values[0])
    if distance_m >= float(distances_m[-1]):
        return float(values[-1])
    return float(np.interp(distance_m, distances_m, values))


def append_line(lines: list[str], label: str, value: str) -> None:
    lines.append(f"{label:<43}: {value}")


def append_selection_process(lines: list[str], payload: dict, config: dict) -> None:
    opt_cfg = optimization_config(config)
    priority_labels = {
        "assist_blower_pressure_excess_Pa": "assist-blower static-pressure excess above the available fan-curve limit (Pa)",
        "current_excess_A": "current excess above the supply limit (A)",
        "wire_temp_excess_C": "wire-temperature excess above the limit (C)",
        "air_outlet_excess_C": "outlet-air excess above the limit (C)",
        "pressure_drop_excess_Pa": "pressure-drop excess above the limit (Pa)",
        "water_loss_excess_kg_day": f"water-loss excess above {float(payload['limits'].get('maxWaterLoss_kg_day', np.nan)):.1f} kg/24 h",
        "bottom_temp_shortfall_C": f"bottom-temperature shortfall below {float(opt_cfg['min_bottom_temp_C']):.1f} C",
        "top_temp_shortfall_C": f"top-temperature shortfall below {float(opt_cfg['min_top_temp_C']):.1f} C",
        "temp_spread_excess_C": f"temperature-spread excess above {float(opt_cfg['spread_limit_C']):.1f} C",
        "heat_shortfall": "heat shortfall, expressed as max(0, dutyNeeded - 1)",
        "minus_min_bed_temp_C": "negative colder-node temperature, so warmer colder-node values rank earlier",
        "minus_mean_bed_temp_C": "negative mean bed temperature, so warmer overall bed values rank earlier",
        "power_W": "total electrical power (W)",
        "temp_spread_C": "raw bottom-top temperature spread (C)",
    }

    lines.append("")
    lines.append("Selection process:")
    lines.append("  One consistent lexicographic rule is used everywhere in this study:")
    lines.append("    1. Vessel-comparison recommended point selection.")
    lines.append("    2. Top candidate ranking for the active configuration.")
    lines.append("    3. Wire-diameter hard-safe selected point at each sampled diameter.")
    lines.append("    4. Optimization-ranked point table.")
    lines.append("  Step 1: compute the hard-constraint residuals listed below for every candidate point.")
    lines.append("  Step 2: define a hard-safe point as one with zero residual for current, wire temperature,")
    lines.append("          outlet-air temperature, pressure drop, water loss, hole velocity, bottom temperature, top temperature,")
    lines.append("          and bottom-top spread, plus positive net heat-to-bed output.")
    lines.append("  Step 3: if one or more hard-safe points exist, rank only those points lexicographically")
    lines.append("          with the ordered criteria below. Inside that hard-safe set, the clean physical")
    lines.append("          rule is: maximize the colder-node bed temperature, then minimize power, then")
    lines.append("          minimize bottom-top spread; mean bed temperature and heat shortfall are only")
    lines.append("          late tie-breakers.")
    lines.append("  Step 4: if no hard-safe point exists, report no recommended design point and fall back to")
    lines.append("          the lexicographically best available point only as a plot/reference point.")
    lines.append("  No blended penalty score is used anywhere in the selection logic.")
    lines.append("  Criteria settings:")
    lines.append(f"    spread limit                       : {float(opt_cfg['spread_limit_C']):.1f} C")
    lines.append(f"    minimum bottom temperature         : {float(opt_cfg['min_bottom_temp_C']):.1f} C")
    lines.append(f"    minimum top temperature            : {float(opt_cfg['min_top_temp_C']):.1f} C")
    lines.append("  Hard-constraint residuals            : current, wire temperature, outlet-air temperature,")
    lines.append("                                         pressure drop, water loss, hole-velocity excess, bottom-temperature shortfall,")
    lines.append("                                         top-temperature shortfall, and spread excess")
    lines.append("  Priority order:")
    for idx, key in enumerate(opt_cfg["priority_order"], start=1):
        lines.append(f"    {idx:>2}. {priority_labels.get(key, key)}")
    lines.append("  Wire-diameter study note:")
    lines.append("    The wire plot now distinguishes hard-safe selected points from best-available")
    lines.append("    reference points. If no hard-safe point exists at a sampled diameter, the solid")
    lines.append("    selected curve is left blank there rather than dropping onto a last-resort point.")


def append_optimization_process(lines: list[str], payload: dict, config: dict) -> None:
    opt_cfg = optimization_config(config)
    if not opt_cfg.get("enabled", True):
        return

    lines.append("")
    lines.append("Optimization process:")
    lines.append("  The optimization table is not a separate penalty-score solver.")
    lines.append("  It is the same lexicographic ranking applied to the active-configuration design-point list.")
    lines.append("  Hard-safe points rank ahead of violating points. If no hard-safe points exist, the table")
    lines.append("  reports the least-violating best-available points and labels the active failure mode.")
    lines.append("  Clean objective inside the hard-safe set:")
    lines.append("    1. maximize min(T_bottom, T_top)")
    lines.append("    2. minimize total electrical power")
    lines.append("    3. minimize |T_bottom - T_top|")
    lines.append("    4. use mean bed temperature and heat shortfall only as late tie-breakers")
    lines.append("  Interpretation for a real system:")
    guidance = opt_cfg.get("real_world_guidance", [])
    if guidance:
        for idx, item in enumerate(guidance, start=1):
            lines.append(f"    {idx}. {item}")
    else:
        lines.append("    1. Keep equipment and safety limits as hard constraints.")
        lines.append("    2. Keep bottom/top temperature minimums and spread limit as process criteria.")
        lines.append("    3. After those are met, prefer warmer bed temperatures and then lower power.")


def append_moisture_process(lines: list[str], payload: dict, config: dict) -> None:
    limits = payload["limits"]
    model_overrides = normalized_model_overrides(config)
    perforation_inputs = model_overrides.get("perforation", {})
    evaporation_inputs = model_overrides.get("evaporation", {})
    cooling_cfg = cooling_mode_config(config)
    aeration_inputs = model_overrides.get("aeration", {})
    calibration_inputs = model_overrides.get("calibration", {})
    summary_inputs = config.get("summary_inputs", {})
    bin_inputs = effective_bin_inputs(summary_inputs, model_overrides)
    wetted = derived_wetted_area_details(summary_inputs, model_overrides)
    wetted_area = float(wetted.get("area_m2", np.nan))
    mode = str(payload.get("meta", {}).get("operatingMode", HEATING_MODE))
    rh_value = float(evaporation_inputs.get("relativeHumidity", np.nan))
    surface_temp_value = float(evaporation_inputs.get("surfaceTemp_C", np.nan))
    if mode == COOLING_MODE:
        rh_value = float(cooling_cfg.get("inletRelativeHumidity", rh_value))
        surface_temp_value = float(cooling_cfg.get("surfaceTemp_C", surface_temp_value))

    lines.append("")
    lines.append("Perforation and moisture process:")
    lines.append("  Per-hole release and moisture loss are evaluated explicitly for the aeration tubes.")
    lines.append("  Source map                           : Frederickson et al. (2007) Table 3 particle-size data,")
    lines.append("                                         Churchill-Bernstein hole convection, Lewis relation,")
    lines.append("                                         and IAPWS saturation/latent-heat relations")
    if {"holesPerTube", "holeDiameter_m"} <= perforation_inputs.keys():
        hole_area_m2 = np.pi * float(perforation_inputs["holeDiameter_m"]) ** 2 / 4.0
        lines.append(
            "  Perforation geometry                 : "
            f"{int(perforation_inputs['holesPerTube'])} holes/tube, "
            f"{1000 * float(perforation_inputs['holeDiameter_m']):.1f} mm diameter, "
            f"area per hole {hole_area_m2:.8f} m^2"
        )
    if "releaseFraction" in aeration_inputs:
        lines.append(
            "  Tube release fraction                : "
            f"{100 * float(aeration_inputs['releaseFraction']):.1f}% of each tube flow is released through the perforations"
        )
    lines.append("  Per-hole flow relation               : released tube flow / holes per tube")
    lines.append("  Hole-velocity relation               : released tube flow / (holes per tube x area per hole)")
    lines.append("  Evaporation relation                 : m_dot = min(m_diffusion, m_capacity, m_heat)")
    lines.append("  Diffusion limit                      : m_diffusion = h_m A (rho_v,sat(T_surface) - rho_v,bulk)")
    lines.append("                                         the diffusion, carrying-capacity, and heat-limited")
    lines.append("                                         evaporation rates are all evaluated in kg/s")
    lines.append("  Heat-mass analogy                    : h_m = h / (rho c_p Le^(2/3)); h comes from the hole-scale convection correlation")
    lines.append("  Sensible-to-latent cap               : m_heat = q_available / h_fg")
    if evaporation_inputs:
        wetted_area_basis = str(evaporation_inputs.get("wettedAreaBasis", "manual")).strip().lower()
        lines.append(
            "  Moisture-state inputs                : "
            f"RH {100 * rh_value:.0f}%, "
            f"surface temp {surface_temp_value:.1f} C, "
            f"wetted area {wetted_area:.3f} m^2, "
            f"Lewis factor {float(evaporation_inputs.get('lewisFactor', np.nan)):.2f}"
        )
        if wetted_area_basis == "top_plan_area":
            wet_fraction = float(evaporation_inputs.get("wettedAreaFraction", 1.0))
            lines.append(
                "  Wetted-area assumption               : A_wet = f_wet x L x W = "
                f"{wet_fraction:.2f} x {float(bin_inputs['length_m']):.2f} x "
                f"{float(bin_inputs['width_m']):.2f} = {wetted_area:.3f} m^2"
            )
        elif wetted_area_basis == "table3_particle_surface":
            labels = [PARTICLE_SURFACE_LABELS[key] for key in PARTICLE_SURFACE_KEYS]
            windrow_pct = wetted.get("windrow_pct", np.full(len(labels), np.nan))
            vermi_pct = wetted.get("vermi_pct", np.full(len(labels), np.nan))
            blended = 100.0 * wetted.get("blended_weights", np.full(len(labels), np.nan))
            reps = wetted.get("representative_mm", np.full(len(labels), np.nan))
            lines.append(
                "  Wetted-area basis                    : particle-surface estimate from Frederickson et al. (2007), Table 3"
            )
            lines.append(
                "  Table 3 windrow fractions            : "
                + ", ".join(f"{label} {value:.1f}%" for label, value in zip(labels, windrow_pct))
            )
            lines.append(
                "  Table 3 vermicompost fractions       : "
                + ", ".join(f"{label} {value:.1f}%" for label, value in zip(labels, vermi_pct))
            )
            lines.append(
                "  Blended fractions used               : "
                + ", ".join(f"{label} {value:.1f}%" for label, value in zip(labels, blended))
            )
            lines.append(
                "  Representative diameters             : "
                + ", ".join(f"{label} {value:.1f} mm" for label, value in zip(labels, reps))
            )
            lines.append(
                "  Equivalent particle diameter         : "
                f"d32 = {1000.0 * float(wetted.get('d32_m', np.nan)):.2f} mm from "
                "1/d32 = sum(w_i / d_i)"
            )
            lines.append(
                "  External area per fill volume        : "
                f"6/d32 = {float(wetted.get('area_per_volume_m2_m3', np.nan)):.1f} m^2/m^3"
            )
            lines.append(
                "  Wetted-area derivation               : "
                f"A_wet = f_access x (6/d32) x V_fill = {float(wetted.get('accessible_fraction', np.nan)):.2f} x "
                f"{float(wetted.get('area_per_volume_m2_m3', np.nan)):.1f} x {float(wetted.get('fill_volume_m3', np.nan)):.3f} = "
                f"{wetted_area:.3f} m^2"
            )
            lines.append("                                         assumes the Table 3 mass fractions represent equal-density")
            lines.append("                                         agglomerate classes, so the external particle area follows")
            lines.append("                                         the Sauter-mean relation for equivalent spheres")
    if "evaporationHMultiplier" in calibration_inputs:
        lines.append(
            "  Evaporation h multiplier             : "
            f"{float(calibration_inputs['evaporationHMultiplier']):.2f}"
        )
    if mode == COOLING_MODE:
        lines.append(
            "  Water-loss constraint                : "
            f"predicted substrate water loss must stay at or below {float(limits.get('maxWaterLoss_kg_day', np.nan)):.1f} kg/24 h"
        )
    else:
        lines.append(
            "  Water-loss constraint                : "
            f"predicted sprinkler make-up must stay at or below {float(limits.get('maxWaterLoss_kg_day', np.nan)):.1f} kg/24 h"
        )
    lines.append("  Velocity-limit interpretation        : there is no single fixed hole-velocity limit independent of")
    lines.append("                                         voltage, temperature, humidity, and flow; the dashed red")
    lines.append("                                         contour on the moisture plot is the actual allowable boundary")
    lines.append("                                         for the configured water-loss limit")
    lines.append("  Bed-side h interpretation            : the configured bed-side h is treated as an effective sensible")
    lines.append("                                         coefficient; latent evaporation is modeled separately and is")
    lines.append("                                         therefore not folded into h a second time")


def append_radial_profile_process(lines: list[str], payload: dict, config: dict) -> None:
    radial = compute_radial_profile(payload, config)
    if not radial:
        return

    lines.append("")
    lines.append("Radial air-cooling profile:")
    lines.append("  Purpose                              : estimate the distance from the aeration tube outer surface")
    lines.append("                                         at which the released air cools to the configured worm-safety")
    lines.append("                                         threshold")
    lines.append("  Source map                           : finite-radius plume energy balance with a local")
    lines.append("                                         temperature-plus-vapor re-solve using the same")
    lines.append("                                         Table 3 wetted-area basis as the moisture model")
    lines.append("  Governing temperature relation       : m'_rel c_p dT/dx = -(q'_sens + q'_lat)")
    lines.append("  Local sensible term                  : q'_sens = h_sens pi a_wet b(x)^2 (T - T_sink)")
    lines.append("  Local evaporation term               : q'_lat = m'_evap h_fg with")
    lines.append("                                         m'_evap = h_m pi a_wet b(x)^2 (rho_v,sat(T_surface) - rho_v,bulk)")
    lines.append("                                         limited stepwise by plume carrying capacity and by the")
    lines.append("                                         remaining sensible enthalpy above T_sink")
    lines.append("  Plume-radius law                     : b(x) = b0 + alpha x, with optional cap b_max")
    lines.append(
        "  Tube outer radius                    : "
        f"{1000.0 * float(radial['tube_outer_radius_m']):.1f} mm"
    )
    lines.append(
        "  Initial plume radius                 : "
        f"{1000.0 * float(radial['initial_plume_radius_m']):.1f} mm"
    )
    lines.append(
        "  Plume spreading rate                 : "
        f"{float(radial['spread_rate_m_per_m']):.3f} m plume radius per m radial distance"
    )
    if np.isfinite(float(radial["max_plume_radius_m"])):
        lines.append(
            "  Maximum plume radius                 : "
            f"{1000.0 * float(radial['max_plume_radius_m']):.1f} mm"
        )
    else:
        lines.append("  Maximum plume radius                 : none; plume radius grows linearly across the plotted range")
    lines.append(f"  Selected outlet temperature          : {float(active_reference_point(payload).get('airOutlet_C', np.nan)):.1f} C")
    lines.append(f"  Sink temperature used                : {float(radial['sink_C']):.1f} C")
    lines.append(f"  Bed reference temperature            : {float(radial['bed_ref_C']):.1f} C")
    lines.append(f"  Wetted-area density                  : {float(radial['area_density_m2_m3']):.1f} m^2/m^3")
    lines.append(f"  Base sensible coefficient            : {float(radial['h_sensible_W_m2K']):.2f} W/m^2-K")
    lines.append(f"  Mean latent coefficient              : {float(radial['h_latent_W_m2K']):.2f} W/m^2-K")
    lines.append(f"  Mean combined sink coefficient       : {float(radial['h_effective_W_m2K']):.2f} W/m^2-K")
    lines.append(
        "  Released mass flow per tube length   : "
        f"{float(radial['mdot_release_per_length_kg_s_m']):.6f} kg/s-m"
    )
    lines.append(
        "  Initial / final plume RH             : "
        f"{100.0 * float(radial['relative_humidity_inlet']):.1f}% -> {100.0 * float(radial['relative_humidity_final']):.1f}%"
    )
    lines.append(
        "  Initial / final vapor density        : "
        f"{float(radial['initial_bulk_vapor_density_kg_m3']):.4f} -> {float(radial['final_bulk_vapor_density_kg_m3']):.4f} kg/m^3"
    )
    lines.append(
        "  Integrated radial sensible / latent  : "
        f"{float(radial['total_radial_sensible_W']):.2f} W and {float(radial['total_radial_latent_W']):.2f} W"
    )
    if radial["threshold_state"] == "not_needed":
        lines.append(
            "  Safe distance result                 : 0.0 mm because the selected outlet temperature is already"
        )
        lines.append(
            f"                                         at or below the {float(radial['threshold_C']):.1f} C threshold"
        )
    elif radial["threshold_state"] == "reached":
        lines.append(
            "  Safe distance result                 : "
            f"{1000.0 * float(radial['safe_distance_m']):.1f} mm from the tube outer surface"
        )
    elif radial["threshold_state"] == "beyond_plot_limit":
        lines.append(
            "  Safe distance result                 : "
            f"{1000.0 * float(radial['safe_distance_m']):.1f} mm; this is beyond the plotted radial range"
        )
    else:
        lines.append("  Safe distance result                 : no finite crossing in this simple model because")
        lines.append(
            f"                                         the profile asymptote stays above {float(radial['threshold_C']):.1f} C"
        )


def append_tube_layout_process(
    lines: list[str],
    payload: dict,
    config: dict,
    layout: dict | None = None,
) -> None:
    if layout is None:
        layout = compute_tube_layout_and_habitat(payload, config)
    if not layout:
        return

    fill_volume = float(layout["fill_volume_m3"])
    exclusion_volume = float(layout["exclusion_volume_m3"])
    pipe_volume = float(layout["pipe_volume_m3"])
    wall_return_volume = float(layout.get("wall_return_volume_m3", 0.0))
    unsafe_volume = float(layout.get("unsafe_volume_m3", exclusion_volume))
    unsafe_shell = float(layout["unsafe_shell_volume_m3"])
    habitable_volume = float(layout["habitable_volume_m3"])
    exclusion_fraction = 100.0 * exclusion_volume / max(fill_volume, 1e-12)
    unsafe_fraction = 100.0 * unsafe_volume / max(fill_volume, 1e-12)
    pipe_fraction = 100.0 * pipe_volume / max(fill_volume, 1e-12)
    habitable_fraction = 100.0 * habitable_volume / max(fill_volume, 1e-12)

    center_mm = ", ".join(f"{1000.0 * value:.1f}" for value in layout["x_centers_m"])
    if layout["pitches_m"].size > 0:
        pitch_mm = ", ".join(f"{1000.0 * value:.1f}" for value in layout["pitches_m"])
        clear_gap_mm = ", ".join(f"{1000.0 * value:.1f}" for value in layout["clear_gaps_m"])
    else:
        pitch_mm = "n/a"
        clear_gap_mm = "n/a"

    lines.append("")
    lines.append("Tube layout and worm-habitat exclusion:")
    if str(layout.get("center_mode", "")).startswith("rule_based"):
        lines.append("  Tube-center placement                : rule-based 1-to-12 tube layout from tube_layout.layoutMode")
        lines.append("                                         using 10% of the filled-bed width/height as the wall offset")
    elif layout.get("center_mode") == "uniform_across_width":
        lines.append("  Tube-center placement                : auto-generated uniform spacing across the bed width")
        lines.append("                                         x_i = (i - 0.5) W / N_tubes")
    elif layout.get("center_mode") == "explicit_xy":
        lines.append("  Tube-center placement                : user-specified (x, y) positions from tube_layout.tubeCenters_m")
    else:
        lines.append("  Tube-center placement                : user-specified x locations from tube_layout.tubeCenterX_m")
    lines.append(f"  Filled-bed cross section             : {1000.0 * layout['width_m']:.1f} mm wide x {1000.0 * layout['fill_height_m']:.1f} mm high")
    lines.append(f"  Tube center height above floor       : {1000.0 * float(layout['y_centers_m'][0]):.1f} mm")
    lines.append(f"  Tube center x positions              : {center_mm} mm")
    lines.append(f"  Center-to-center pitch               : {pitch_mm} mm")
    lines.append(f"  Clear gap between tube OD surfaces   : {clear_gap_mm} mm")
    if np.isfinite(float(layout.get("nearest_center_spacing_m", np.nan))):
        lines.append(
            "  Nearest tube-center spacing          : "
            f"{1000.0 * float(layout['nearest_center_spacing_m']):.1f} mm"
        )
        lines.append(
            "  Nearest clear tube-to-tube gap       : "
            f"{1000.0 * float(layout['nearest_clear_gap_m']):.1f} mm"
        )
    lines.append(
        "  Side clearances to tube OD           : "
        f"left {1000.0 * layout['side_clearance_left_m']:.1f} mm, right {1000.0 * layout['side_clearance_right_m']:.1f} mm"
    )
    lines.append(
        "  Safe-distance exclusion radius       : "
        f"tube radius + safe offset = {1000.0 * layout['tube_radius_m']:.1f} + {1000.0 * layout['safe_distance_m']:.1f} = "
        f"{1000.0 * layout['exclusion_radius_m']:.1f} mm"
    )
    lines.append(
        "  Pipe solid volume                    : "
        f"{pipe_volume:.4f} m^3 ({pipe_fraction:.2f}% of filled-bed volume)"
    )
    lines.append(
        "  Pipe exclusion-cylinder union volume : "
        f"{exclusion_volume:.4f} m^3 ({exclusion_fraction:.2f}% of filled-bed volume)"
    )
    lines.append(
        "  Added wall-return unsafe volume      : "
        f"{wall_return_volume:.4f} m^3 before overlap with the tube exclusion cylinders"
    )
    lines.append(
        "  Total worm-unsafe union volume       : "
        f"{unsafe_volume:.4f} m^3 ({unsafe_fraction:.2f}% of filled-bed volume)"
    )
    lines.append(
        "  Worm-unsafe shell volume             : "
        f"{unsafe_shell:.4f} m^3 ({100.0 * unsafe_shell / max(fill_volume, 1e-12):.2f}% of filled-bed volume)"
    )
    lines.append(
        "  Remaining habitable bed volume       : "
        f"{habitable_volume:.4f} m^3 ({habitable_fraction:.2f}% of filled-bed volume)"
    )
    lines.append(
        "  Pipe-overlap correction              : "
        f"{float(layout['pipe_overlap_volume_m3']):.4f} m^3 removed from the naive sum of {int(layout['n_tubes'])} pipe cylinders"
    )
    lines.append(
        "  Exclusion-overlap correction         : "
        f"{float(layout['exclusion_overlap_volume_m3']):.4f} m^3 removed from the naive sum of {int(layout['n_tubes'])}"
    )
    lines.append("                                         safe-distance exclusion cylinders")
    wall_return = dict(layout.get("wall_return", {}))
    if any(float(wall_return.get(f"{side}_depth_m", 0.0)) > 0 for side in ("left", "right", "bottom")):
        lines.append(
            "  Wall-return unsafe depths            : "
            f"left {1000.0 * float(wall_return.get('left_depth_m', 0.0)):.1f} mm, "
            f"right {1000.0 * float(wall_return.get('right_depth_m', 0.0)):.1f} mm, "
            f"bottom {1000.0 * float(wall_return.get('bottom_depth_m', 0.0)):.1f} mm"
        )
        lines.append(
            "  Wall plume-contact temperatures      : "
            f"left {float(wall_return.get('left_contactTemp_C', np.nan)):.1f} C, "
            f"right {float(wall_return.get('right_contactTemp_C', np.nan)):.1f} C, "
            f"bottom {float(wall_return.get('bottom_contactTemp_C', np.nan)):.1f} C"
        )
        lines.append(
            "  Predicted inner-wall temperatures    : "
            f"left {float(wall_return.get('left_wallTemp_C', np.nan)):.1f} C, "
            f"right {float(wall_return.get('right_wallTemp_C', np.nan)):.1f} C, "
            f"bottom {float(wall_return.get('bottom_wallTemp_C', np.nan)):.1f} C"
        )
    lines.append(
        "  Overlap evaluation method            : overlap-aware union area on a "
        f"{int(layout['cross_section_grid_points'])} x {int(layout['cross_section_grid_points'])} midpoint grid"
    )
    lines.append("                                         clipped to the filled-bed rectangle, then extruded over the")
    lines.append("                                         aeration-tube length to obtain volume; the final unsafe union")
    lines.append("                                         includes both plume exclusion cylinders and 1D wall-return bands")


def normalize_label(label: str) -> str:
    return label.strip().lower().replace(" ", "_").replace("+", "plus").replace("-", "_")


def display_configuration_label(label: str) -> str:
    mapping = {
        "covered_top_plus_vent": "closed top with vents",
        "covered_top": "closed top",
        "uncovered_top": "open top",
        "open_top": "open top",
    }
    normalized = normalize_label(label)
    return mapping.get(normalized, label.strip().lower())


def same_configuration_label(active_label: str, study_label: str) -> bool:
    a = normalize_label(active_label)
    b = normalize_label(study_label)
    if a == b:
        return True
    aliases = {
        ("open_top", "uncovered_top"),
    }
    return (a, b) in aliases or (b, a) in aliases


def configuration_extras(config: dict, label: str) -> dict:
    return config.get("summary_inputs", {}).get("configuration_reports", {}).get(normalize_label(label), {})


def trade_entry_valid_mask(entry: dict, config: dict) -> np.ndarray:
    x_values = as_float_array(entry.get("xValues", []))
    if x_values.size == 0:
        return np.array([], dtype=bool)

    mask = np.isfinite(x_values)
    if str(entry.get("xField", "")) == "coilMeanD_mm":
        model_overrides = normalized_model_overrides(config)
        heater_inputs = effective_heater_inputs(config.get("summary_inputs", {}), model_overrides)
        wire_inputs = config.get("summary_inputs", {}).get("wire", {})
        heater_id_mm = float(heater_inputs.get("ID_mm", np.nan))
        wire_diameter_mm = float(wire_inputs.get("diameter_mm", np.nan))
        max_fit_mm = heater_id_mm - wire_diameter_mm
        if np.isfinite(max_fit_mm):
            mask &= x_values <= max_fit_mm + 1e-9
    return mask


def choose_trade_display_index(entry: dict, config: dict, selected: bool) -> int | None:
    x_values = as_float_array(entry.get("xValues", []))
    if x_values.size == 0:
        return None

    valid_mask = trade_entry_valid_mask(entry, config)
    if not np.any(valid_mask):
        return None

    index_key = "bestIndex" if selected else "bestAvailableIndex"
    raw_best = entry.get(index_key)
    if raw_best is not None and np.isfinite(float(raw_best)):
        candidate = int(float(raw_best)) - 1
        if 0 <= candidate < len(x_values) and valid_mask[candidate]:
            return candidate

    if selected:
        primary = as_float_array(entry.get("selectedColderNode_C", []))
        secondary = as_float_array(entry.get("selectedPower_W", []))
        tertiary = as_float_array(entry.get("selectedCurrent_A", []))
    else:
        primary = as_float_array(entry.get("bestAvailableColderNode_C", []))
        secondary = as_float_array(entry.get("bestAvailablePower_W", []))
        tertiary = as_float_array(entry.get("bestAvailableCurrent_A", []))

    if primary.size != len(x_values):
        primary = np.full(len(x_values), np.nan)
    if secondary.size != len(x_values):
        secondary = np.full(len(x_values), np.nan)
    if tertiary.size != len(x_values):
        tertiary = np.full(len(x_values), np.nan)

    valid_idx = np.where(valid_mask & np.isfinite(primary) & np.isfinite(secondary) & np.isfinite(tertiary))[0]
    if valid_idx.size == 0 and selected:
        return None
    if valid_idx.size == 0:
        valid_idx = np.where(valid_mask)[0]
    if valid_idx.size == 0:
        return None

    order = np.lexsort((tertiary[valid_idx], secondary[valid_idx], -primary[valid_idx]))
    return int(valid_idx[order[0]])


def append_geometry_trade_summary(lines: list[str], payload: dict, config: dict) -> None:
    studies = [
        ("wire_diameter_sweep", "Wire diameter study"),
        ("coil_diameter_sweep", "Coil mean diameter study"),
        ("wire_length_sweep", "Wire length study"),
    ]
    active_label = str(payload.get("active_configuration", {}).get("label", ""))
    study_lines: list[str] = []
    for key, title in studies:
        entries = payload.get(key, [])
        entry = next(
            (item for item in entries if same_configuration_label(active_label, item.get("label", ""))),
            None,
        )
        if not entry:
            continue
        x_values = as_float_array(entry.get("xValues", []))
        selected_idx = choose_trade_display_index(entry, config, selected=True)
        available_idx = choose_trade_display_index(entry, config, selected=False)
        if selected_idx is not None and 0 <= selected_idx < len(x_values):
            study_lines.append(
                f"  {title:<36}: hard-safe optimum at {x_values[selected_idx]:.3f} {entry.get('xUnit', '')}, "
                f"colder node {float(as_float_array(entry.get('selectedColderNode_C', []))[selected_idx]):.1f} C, "
                f"Qbed {float(as_float_array(entry.get('selectedQtoBed_W', []))[selected_idx]):.1f} W"
            )
        elif available_idx is not None and 0 <= available_idx < len(x_values):
            study_lines.append(
                f"  {title:<36}: best-available only at {x_values[available_idx]:.3f} {entry.get('xUnit', '')}, "
                f"colder node {float(as_float_array(entry.get('bestAvailableColderNode_C', []))[available_idx]):.1f} C, "
                f"Qbed {float(as_float_array(entry.get('bestAvailableQtoBed_W', []))[available_idx]):.1f} W"
            )
    if study_lines:
        lines.append("")
        lines.append("Helical geometry trade studies:")
        lines.extend(study_lines)


def write_design_point_block(
    lines: list[str],
    heading: str,
    config_label: str,
    rec: dict,
    summary_inputs: dict,
    limits: dict,
    opt_cfg: dict,
    sweep: dict,
    extras: dict,
) -> None:
    heater_inputs = summary_inputs.get("heater_tube", {})
    aeration_inputs = summary_inputs.get("aeration", {})
    wire_inputs = summary_inputs.get("wire", {})

    lines.append("")
    lines.append(heading)
    lines.append(f"  Configuration                      : {display_configuration_label(config_label)}")
    if not rec:
        lines.append("  Result                             : no hard-constraint-safe design point found in the current sweep")
        return

    n_tubes = int(aeration_inputs.get("nParallelTubes", 1))
    branch_flow_lpm = float(rec["totalFlow_Lpm"]) / max(n_tubes, 1)
    branch_current_a = float(rec["totalCurrent_A"]) / max(n_tubes, 1)
    branch_power_w = float(rec["totalPower_W"]) / max(n_tubes, 1)
    branch_resistance_ohm = float(rec["voltage_V"]) / max(branch_current_a, 1e-9)
    supply_eq_resistance_ohm = float(rec["voltage_V"]) / max(float(rec["totalCurrent_A"]), 1e-9)
    wire_length_m = rec.get("wireLength_m")
    if wire_length_m is None or not np.isfinite(float(wire_length_m)):
        wire_length_m = compute_wire_length_m(float(rec["pitch_mm"]), heater_inputs, wire_inputs)

    lines.append(f"  Pitch                              : {float(rec['pitch_mm']):.1f} mm")
    lines.append(f"  Supply voltage per branch          : {float(rec['voltage_V']):.1f} V")
    lines.append(f"  Total flow                         : {float(rec['totalFlow_Lpm']):.1f} L/min")
    lines.append(
        "  Selected operating trio            : "
        f"{float(rec['totalPower_W']):.1f} W, {float(rec['totalCurrent_A']):.2f} A, {float(rec['totalFlow_Lpm']):.1f} L/min"
    )
    lines.append(f"  Branch flow                        : {branch_flow_lpm:.1f} L/min")
    if wire_length_m is not None:
        lines.append(f"  Branch wire length                 : {wire_length_m:.2f} m")
    if "coilSpan_m" in heater_inputs:
        lines.append(f"  Branch coil axial length           : {float(heater_inputs['coilSpan_m']):.2f} m")
    if {"gauge", "diameter_mm"} <= wire_inputs.keys():
        lines.append(f"  Wire gauge / diameter              : {wire_inputs['gauge']} / {float(wire_inputs['diameter_mm']):.3f} mm")
    lines.append(f"  Branch resistance                  : {branch_resistance_ohm:.2f} ohm")
    lines.append(f"  Supply equivalent resistance       : {supply_eq_resistance_ohm:.2f} ohm")
    lines.append(f"  Branch current                     : {branch_current_a:.2f} A")
    lines.append(f"  Total current                      : {float(rec['totalCurrent_A']):.2f} A")
    lines.append(f"  Branch electrical power            : {branch_power_w:.1f} W")
    lines.append(f"  Total electrical power             : {float(rec['totalPower_W']):.1f} W")
    if rec.get("hWire_W_m2K") is not None:
        lines.append(f"  h_wire                             : {float(rec['hWire_W_m2K']):.1f} W/m^2-K")
    elif extras.get("hWire_W_m2K") is not None:
        lines.append(f"  h_wire                             : {float(extras['hWire_W_m2K']):.1f} W/m^2-K")
    else:
        lines.append("  h_wire                             : not provided in exported JSON/config")
    lines.append(f"  Branch air outlet temperature      : {float(rec['airOutlet_C']):.1f} C")
    lines.append(f"  Max wire temperature               : {float(rec['wireMax_C']):.1f} C")
    lines.append(f"  Total heat to bed                  : {float(rec['QtoBed_W']):.1f} W")
    lines.append(f"  Duty needed at design loss         : {float(rec['dutyNeeded']):.2f}")
    lines.append(f"  Total pressure drop                : {float(rec['deltaP_Pa']):.1f} Pa")
    if np.isfinite(float(rec.get("distributionDeltaP_Pa", np.nan))):
        lines.append(f"  Distribution manifold dP          : {float(rec['distributionDeltaP_Pa']):.1f} Pa")
        lines.append(
            "  Distribution dP components         : "
            f"header {float(rec.get('distributionHeader_Pa', 0.0)):.1f}, "
            f"header->splitter {float(rec.get('distributionHeaderToSplitter_Pa', 0.0)):.1f}, "
            f"splitter body {float(rec.get('distributionSplitterBody_Pa', 0.0)):.1f}, "
            f"connector friction {float(rec.get('distributionConnectorFriction_Pa', 0.0)):.1f}, "
            f"connector->branch {float(rec.get('distributionConnectorToBranch_Pa', 0.0)):.1f} Pa"
        )
    if "flowPerHole_Lpm" in rec:
        lines.append(f"  Released flow per hole             : {float(rec['flowPerHole_Lpm']):.2f} L/min")
    if "holeVelocity_m_s" in rec:
        lines.append(f"  Hole velocity                      : {float(rec['holeVelocity_m_s']):.2f} m/s")
    if "waterLoss_kg_day" in rec:
        lines.append(f"  Water loss / sprinkler make-up     : {float(rec['waterLoss_kg_day']):.2f} kg/24 h")
    if "latentEvap_W" in rec:
        lines.append(f"  Latent evaporation load            : {float(rec['latentEvap_W']):.1f} W")
    bottom_full = rec.get("TbottomFullPower_C", extras.get("fullPowerBottom_C"))
    top_full = rec.get("TtopFullPower_C", extras.get("fullPowerTop_C"))
    if bottom_full is not None and top_full is not None:
        lines.append(
            "  Full-power bed estimate            : "
            f"bottom {float(bottom_full):.1f} C, top {float(top_full):.1f} C"
        )
    else:
        lines.append("  Full-power bed estimate            : not provided in exported JSON/config")
    lines.append(f"  Status                             : {recommended_status_long(rec, limits, opt_cfg)}")
    if float(rec["dutyNeeded"]) > float(sweep["targetDuty"]):
        lines.append(f"  Note                               : {100 * float(sweep['targetDuty']):.0f}% duty is not achieved;")
        lines.append("                                       more tubes help local temperature")
        lines.append("                                       distribution, but not the minimum")
        lines.append("                                       m_dot*cp*DeltaT requirement")


def cooling_status_long(rec: dict) -> str:
    if not rec:
        return "no hard-constraint-safe summer cooling point exists in the current sweep"
    if bool(rec.get("isConstraintSafeByCriteria", False)):
        return "meets all summer hard constraints and provides net bed cooling"
    return "best-available only; it violates one or more summer hard constraints"


def write_cooling_design_point_block(lines: list[str], heading: str, rec: dict, config: dict) -> None:
    lines.append("")
    lines.append(heading)
    if not rec:
        lines.append("  Result                             : no hard-constraint-safe summer cooling point found in the current sweep")
        return

    limits = cooling_mode_config(config)["limits"]
    lines.append(f"  Total flow                         : {float(rec['totalFlow_Lpm']):.1f} L/min")
    lines.append(f"  Spot-cooler target supply air      : {float(rec.get('spotCoolerTargetSupply_C', np.nan)):.1f} C")
    if np.isfinite(float(rec.get("spotCoolerDerivedTargetSupply_C", np.nan))):
        lines.append(f"  Spec-derived middle-case supply    : {float(rec['spotCoolerDerivedTargetSupply_C']):.1f} C")
    lines.append(f"  Spot-cooler actual supply air      : {float(rec.get('spotCoolerActualSupply_C', rec.get('airInlet_C', np.nan))):.1f} C")
    lines.append(f"  Spot-cooler supply RH              : {100.0 * float(rec.get('spotCoolerSupplyRelativeHumidity', np.nan)):.1f} %")
    lines.append(f"  Spot-cooler sensible load          : {float(rec.get('spotCoolerLoad_W', np.nan)):.1f} W")
    if np.isfinite(float(rec.get("spotCoolerPower_W", np.nan))):
        lines.append(f"  Spot-cooler electric power         : {float(rec['spotCoolerPower_W']):.1f} W")
    if np.isfinite(float(rec.get("assistBlowerPower_W", np.nan))):
        lines.append(f"  Assist-blower electric power       : {float(rec['assistBlowerPower_W']):.1f} W")
    if np.isfinite(float(rec.get("assistBlowerAvailablePressure_Pa", np.nan))):
        lines.append(f"  Assist-blower available static     : {float(rec['assistBlowerAvailablePressure_Pa']):.1f} Pa")
    if np.isfinite(float(rec.get("spotCoolerCondensate_kg_day", np.nan))):
        lines.append(f"  Spot-cooler condensate output      : {float(rec['spotCoolerCondensate_kg_day']):.2f} kg/24 h")
    lines.append(f"  Aeration outlet temperature        : {float(rec['airOutlet_C']):.1f} C")
    lines.append(f"  Net cooling to bed                 : {max(-float(rec['QtoBed_W']), 0.0):.1f} W")
    lines.append(f"  Signed Q_to_bed                    : {float(rec['QtoBed_W']):.1f} W")
    lines.append(f"  Bottom-node steady temperature     : {float(rec['TbottomFullPower_C']):.1f} C")
    lines.append(f"  Top-node steady temperature        : {float(rec['TtopFullPower_C']):.1f} C")
    lines.append(f"  Bottom-top spread                  : {abs(float(rec['TbottomFullPower_C']) - float(rec['TtopFullPower_C'])):.1f} C")
    lines.append(f"  Released flow per hole             : {float(rec['flowPerHole_Lpm']):.2f} L/min")
    lines.append(f"  Hole velocity                      : {float(rec['holeVelocity_m_s']):.2f} m/s")
    lines.append(f"  Substrate water loss               : {float(rec['waterLoss_kg_day']):.2f} kg/24 h")
    lines.append(f"  Latent evaporation load            : {float(rec['latentEvap_W']):.1f} W")
    lines.append(f"  Total pressure drop                : {float(rec['deltaP_Pa']):.1f} Pa")
    if np.isfinite(float(rec.get("distributionDeltaP_Pa", np.nan))):
        lines.append(f"  Distribution manifold dP          : {float(rec['distributionDeltaP_Pa']):.1f} Pa")
        lines.append(
            "  Distribution dP components         : "
            f"header {float(rec.get('distributionHeader_Pa', 0.0)):.1f}, "
            f"header->splitter {float(rec.get('distributionHeaderToSplitter_Pa', 0.0)):.1f}, "
            f"splitter body {float(rec.get('distributionSplitterBody_Pa', 0.0)):.1f}, "
            f"connector friction {float(rec.get('distributionConnectorFriction_Pa', 0.0)):.1f}, "
            f"connector->branch {float(rec.get('distributionConnectorToBranch_Pa', 0.0)):.1f} Pa"
        )
    if np.isfinite(float(rec.get("assistBlowerAvailablePressure_Pa", np.nan))):
        lines.append(
            f"  Assist-blower static-pressure margin: "
            f"{float(rec['assistBlowerAvailablePressure_Pa']) - float(rec['deltaP_Pa']):.1f} Pa"
        )
    lines.append(f"  Status                             : {cooling_status_long(rec)}")


def write_cooling_summary(payload: dict, config: dict, output_dir: Path) -> None:
    summary_inputs = config.get("summary_inputs", {})
    model_overrides = normalized_model_overrides(config)
    cooling_cfg = cooling_mode_config(config)
    env_inputs = effective_environment_inputs(summary_inputs, model_overrides)
    env_inputs["greenhouseAir_C"] = float(cooling_cfg["ambientAir_C"])
    env_inputs["pressure_Pa"] = float(cooling_cfg.get("pressure_Pa", env_inputs.get("pressure_Pa", 101325.0)))
    bin_inputs = effective_bin_inputs(summary_inputs, model_overrides)
    wall_inputs = effective_bin_wall_inputs(summary_inputs, model_overrides)
    heater_inputs = effective_heater_inputs(summary_inputs, model_overrides)
    aeration_inputs = effective_aeration_inputs(summary_inputs, model_overrides)
    perforation_inputs = model_overrides.get("perforation", {})
    evaporation_inputs = model_overrides.get("evaporation", {})
    tube_layout = compute_tube_layout_and_habitat(payload, config)
    wetted_area = derived_wetted_area_m2(summary_inputs, model_overrides)
    rec = payload.get("recommended", {})
    best_available = payload.get("best_available", {})
    ranked = rank_cooling_points(payload, config)
    hard_safe_count = sum(1 for pt in ranked if bool(pt.get("isConstraintSafeByCriteria", False)))
    top_list = cooling_top_candidates(payload, config)
    optimization_list = cooling_optimization_ranked_points(payload, config)
    opt = cooling_cfg["optimization"]
    limits = cooling_cfg["limits"]
    lumped_loss = dict(payload.get("lumped_loss", {}))

    lines: list[str] = []
    lines.append("============================================================")
    lines.append("VERMICOMPOSTER SUMMER SPOT-COOLER + ASSIST-BLOWER STUDY")
    lines.append("============================================================")
    append_line(lines, "Operating mode", "summer plenum-fed spot cooler + assist blower")
    append_line(lines, "Hot-weather ambient air temperature", f"{float(cooling_cfg['ambientAir_C']):.1f} C")
    append_line(lines, "Active vessel configuration", "closed top with vents")
    lines.append("")
    lines.append("Design specifications:")
    lines.append(
        "  Vessel interior                         : "
        f"{float(bin_inputs['length_m']):.2f} m x {float(bin_inputs['width_m']):.2f} m x {float(bin_inputs['height_m']):.2f} m, "
        f"fill fraction {float(bin_inputs.get('fillFraction', 0.80)):.2f}"
    )
    vent_count = int(bin_inputs.get("ventHoleCount", 0))
    vent_diameter_m = float(bin_inputs.get("ventHoleDiameter_m", 0.0))
    if vent_count > 0 and vent_diameter_m > 0:
        vent_area = np.pi * vent_diameter_m**2 / 4.0
        lines.append(
            "  Vent specification                     : "
            f"{vent_count} vent hole(s), {1000.0 * vent_diameter_m:.1f} mm diameter each, "
            f"area per hole {vent_area:.6f} m^2, total vent area {vent_count * vent_area:.6f} m^2"
        )
    lines.append(
        "  Vessel wall stack                       : "
        f"{float(wall_inputs['sheetThickness_mm']):.1f} mm steel (k = {float(wall_inputs['sheetK_W_mK']):.1f} W/m-K) + "
        f"{float(wall_inputs['insulationThickness_mm']):.1f} mm insulation (k = {float(wall_inputs['insulationK_W_mK']):.3f} W/m-K)"
    )
    lines.append(
        "  Heater tube hardware                    : "
        f"L {float(heater_inputs['length_m']):.3f} m, ID {float(heater_inputs['ID_mm']):.1f} mm, "
        f"OD {float(heater_inputs['OD_mm']):.1f} mm, insulation {float(heater_inputs['insulationThickness_mm']):.1f} mm"
    )
    lines.append(
        "  Aeration network                        : "
        f"{int(aeration_inputs['nParallelTubes'])} tubes, each L {float(aeration_inputs['length_m']):.2f} m, "
        f"ID {float(aeration_inputs['ID_mm']):.1f} mm, OD {float(aeration_inputs['OD_mm']):.1f} mm"
    )
    header_enabled = bool(aeration_inputs.get("headerEnabled", True))
    splitter_inlet_mm = float(aeration_inputs.get("splitterInletID_mm", aeration_inputs.get("branchConnectorID_mm", 19.0)))
    connector_id_mm = float(aeration_inputs.get("branchConnectorID_mm", 19.0))
    connector_length_m = float(aeration_inputs.get("branchConnectorLength_m", 0.0))
    splitter_outlets = int(round(float(aeration_inputs.get("splitterOutletCount", 4))))
    if header_enabled:
        lines.append(
            "  Distribution manifold                  : "
            f"header ID {float(aeration_inputs.get('headerID_mm', np.nan)):.1f} mm feeding "
            f"{splitter_outlets}-way splitter(s), splitter inlet ID {splitter_inlet_mm:.1f} mm, "
            f"branch connector ID {connector_id_mm:.1f} mm, connector length {connector_length_m:.3f} m"
        )
    else:
        lines.append(
            "  Distribution manifold                  : "
            f"no separate header; {splitter_outlets}-way splitter(s), splitter inlet ID {splitter_inlet_mm:.1f} mm, "
            f"branch connector ID {connector_id_mm:.1f} mm, connector length {connector_length_m:.3f} m"
        )
    lines.append(
        "  Contraction model                      : "
        f"Idel'chik conical converging bellmouth without end wall, "
        f"{float(aeration_inputs.get('contractionConeAngle_deg', 60.0)):.0f} deg"
    )
    if tube_layout:
        lines.append(
            "  Tube layout                             : "
            f"nearest center spacing {1000.0 * float(tube_layout.get('nearest_center_spacing_m', np.nan)):.1f} mm, "
            f"nearest clear OD gap {1000.0 * float(tube_layout.get('nearest_clear_gap_m', np.nan)):.1f} mm"
        )
    if {"holesPerTube", "holeDiameter_m"} <= perforation_inputs.keys():
        hole_area = np.pi * float(perforation_inputs["holeDiameter_m"]) ** 2 / 4.0
        total_holes = int(aeration_inputs.get("nParallelTubes", 1)) * int(perforation_inputs["holesPerTube"])
        lines.append(
            "  Perforations                            : "
            f"{int(perforation_inputs['holesPerTube'])} holes/tube, {1000 * float(perforation_inputs['holeDiameter_m']):.1f} mm dia, "
            f"total holes {total_holes}, total open area {total_holes * hole_area:.6f} m^2"
        )
    spot_cfg = dict(cooling_cfg.get("spot_cooler", {}))
    max_capacity = spot_cfg.get("maxCoolingCapacity_W", None)
    cop = spot_cfg.get("COP", None)
    max_capacity_txt = "unlimited by config" if max_capacity in (None, "") else f"{float(max_capacity):.1f} W"
    cop_txt = "not specified" if cop in (None, "") else f"{float(cop):.2f}"
    target_supply_raw = spot_cfg.get("targetSupplyAir_C", None)
    target_supply_txt = "derived from rated specs" if target_supply_raw in (None, "") else f"{float(target_supply_raw):.1f} C"
    rated_total_btu = spot_cfg.get("ratedTotalCapacity_BTU_h", None)
    rated_airflow = spot_cfg.get("ratedAirflow_CFM", None)
    rated_dehum = spot_cfg.get("ratedDehumidification_pints_day", None)
    derive_from_specs = bool(spot_cfg.get("deriveTargetFromRatedSpecs", False))
    spot_ref = ""
    if best_available:
        spot_ref = str(best_available.get("spotCoolerReferenceNote", ""))
    lines.append(
        "  Spot cooler                            : "
        f"target supply air {target_supply_txt}, "
        f"max sensible capacity {max_capacity_txt}, COP {cop_txt}"
    )
    assist_cfg = dict(cooling_cfg.get("assist_blower", {}))
    lines.append(
        "  Cooling-air architecture               : "
        "spot cooler discharges into an insulated metal plenum; a downstream assist blower draws from the plenum and pushes air through the bed network"
    )
    lines.append(
        "  Assist blower                          : "
        f"{assist_cfg.get('model', 'not specified')}, "
        f"{float(assist_cfg.get('ratedFlow_CFM', np.nan)):.1f} CFM @ "
        f"{float(assist_cfg.get('ratedPressure_inH2O', np.nan)):.1f} in. H2O, "
        f"shutoff {float(assist_cfg.get('shutoffPressure_inH2O', np.nan)):.1f} in. H2O, "
        f"{float(assist_cfg.get('motorPower_W', np.nan)):.1f} W"
    )
    if rated_total_btu not in (None, "") or rated_airflow not in (None, "") or rated_dehum not in (None, ""):
        lines.append(
            "  Spot-cooler rating basis               : "
            f"{float(rated_total_btu) if rated_total_btu not in (None, '') else float('nan'):.0f} BTU/h total, "
            f"{float(rated_airflow) if rated_airflow not in (None, '') else float('nan'):.0f} CFM rated airflow, "
            f"{float(rated_dehum) if rated_dehum not in (None, '') else float('nan'):.1f} pints/day dehumidification"
        )
    if derive_from_specs and best_available:
        lines.append(
            "  Practical middle-case supply estimate  : "
            f"{float(best_available.get('spotCoolerDerivedTargetSupply_C', np.nan)):.1f} C from "
            f"{float(best_available.get('spotCoolerRatedSensibleCapacity_W', np.nan)):.0f} W sensible / "
            f"{float(best_available.get('spotCoolerRatedAirflow_CFM', np.nan)):.0f} CFM at "
            f"{float(cooling_cfg['ambientAir_C']):.1f} C ambient"
        )
        lines.append(
            "  Rated sensible-latent split            : "
            f"{float(best_available.get('spotCoolerRatedSensibleCapacity_W', np.nan)):.0f} W sensible, "
            f"{float(best_available.get('spotCoolerRatedLatentCapacity_W', np.nan)):.0f} W latent"
        )
        lines.append(
            "  Ambient dew-point check                : "
            f"{float(best_available.get('spotCoolerAmbientDewPoint_C', np.nan)):.1f} C dew point at "
            f"{100.0 * float(cooling_cfg.get('inletRelativeHumidity', 0.0)):.0f}% RH; "
            f"{'no condensate expected at the derived supply temperature' if float(best_available.get('spotCoolerCondensate_kg_day', np.nan)) < 1.0e-6 else 'condensate expected during cooling'}"
        )
    if spot_ref:
        lines.append(f"  Spot-cooler source note               : {spot_ref}")
    assist_ref = str(best_available.get("assistBlowerReferenceNote", "")) if best_available else ""
    if assist_ref:
        lines.append(f"  Assist-blower source note             : {assist_ref}")
    lines.append(
        "  Bed moisture model                      : "
        f"ambient RH {100 * float(cooling_cfg.get('inletRelativeHumidity', evaporation_inputs.get('relativeHumidity', 0.0))):.0f}%, "
        f"surface temp {float(cooling_cfg.get('surfaceTemp_C', evaporation_inputs.get('surfaceTemp_C', np.nan))):.1f} C, "
        f"wetted area {wetted_area:.3f} m^2"
    )
    lines.append(
        "  Cooling sweep                           : "
        f"flow = {format_sweep(payload['sweep']['totalFlow_Lpm'])} L/min, electrical heater power fixed at 0 W"
    )
    lines.append(
        "  Summer hard constraints                 : "
        "assist-blower available static >= network dP at the selected flow, "
        f"dP <= {float(limits['maxPressureDrop_Pa']):.0f} Pa, "
        f"water loss <= {float(limits['maxWaterLoss_kg_day']):.1f} kg/24 h, "
        f"hole velocity <= {float(limits['maxHoleVelocity_m_s']):.2f} m/s, "
        f"T_bottom <= {float(opt['max_bottom_temp_C']):.1f} C, "
        f"T_top <= {float(opt['max_top_temp_C']):.1f} C, "
        f"|T_bottom - T_top| <= {float(opt['spread_limit_C']):.1f} C"
    )
    append_line(lines, f"Lumped cooling requirement at {float(cooling_cfg['designBed_C']):.1f} C bed", f"{float(payload['requiredCooling_W']):.1f} W")
    append_line(lines, "Overall UA to ambient air", f"{float(lumped_loss.get('UA_W_K', np.nan)):.3f} W/K")
    append_line(lines, "External natural-convection h estimate", f"{float(lumped_loss.get('externalH_W_m2K', np.nan)):.2f} W/m^2-K")
    append_line(lines, "Lumped time constant", f"{float(lumped_loss.get('tau_h', np.nan)):.1f} h")
    append_line(lines, "Lumped Biot number estimate", f"{float(lumped_loss.get('Bi_lumped', np.nan)):.2f}")
    append_line(lines, "Hard-safe cooling points found", f"{hard_safe_count} of {len(payload['design_points'])}")

    lines.append("")
    lines.append("Summer cooling selection process:")
    lines.append("  Hard-safe means all of the following hold simultaneously:")
    lines.append(
        "    1. the predicted network pressure drop stays below the assist-blower available static pressure at that airflow;"
    )
    lines.append(
        "    2. pressure drop, substrate water loss, and hole velocity remain below their configured limits;"
    )
    lines.append(
        "    3. bottom and top bed temperatures remain below their configured summer maxima;"
    )
    lines.append("    4. bottom-top spread remains below the configured spread limit; and")
    lines.append("    5. net Q_to_bed is negative, so the spot-cooled airflow is actually removing heat from the bed.")
    lines.append("  Clean ordering inside the hard-safe set:")
    lines.append("    1. minimize max(T_bottom, T_top)")
    lines.append("    2. minimize mean bed temperature")
    lines.append("    3. minimize total cooling-system electric power")
    lines.append("    4. minimize total airflow")
    lines.append("    5. minimize bottom-top spread")
    lines.append("  No weighted penalty score is used.")

    append_moisture_process(lines, payload, config)
    append_radial_profile_process(lines, payload, config)
    append_tube_layout_process(lines, payload, config, tube_layout)

    write_cooling_design_point_block(lines, "Criteria-based summer cooling point:", rec, config)
    write_cooling_design_point_block(lines, "Best-available summer reference point:", best_available, config)

    lines.append("")
    if rec:
        lines.append("Top summer cooling candidates:")
    else:
        lines.append("Top best-available summer cooling candidates:")
    lines.append("  Qtot(L/min)  Tbottom(C)  Ttop(C)  dT(C)  Tair,out(C)  cooling(W)  water(kg/d)  Qspot(W)  Pspot(W)  Pblow(W)  Vhole  dP(Pa)  state")
    for pt in top_list:
        spot_power = float(pt.get("spotCoolerPower_W", np.nan))
        blower_power = float(pt.get("assistBlowerPower_W", np.nan))
        lines.append(
            f"  {float(pt['totalFlow_Lpm']):10.1f} {float(pt['TbottomFullPower_C']):10.1f} {float(pt['TtopFullPower_C']):8.1f} "
            f"{abs(float(pt['TbottomFullPower_C']) - float(pt['TtopFullPower_C'])):6.1f} {float(pt['airOutlet_C']):12.1f} "
            f"{max(-float(pt['QtoBed_W']), 0.0):10.1f} {float(pt['waterLoss_kg_day']):11.2f} {float(pt.get('spotCoolerLoad_W', np.nan)):9.1f} "
            f"{spot_power if np.isfinite(spot_power) else float('nan'):9.1f} "
            f"{blower_power if np.isfinite(blower_power) else float('nan'):9.1f} "
            f"{float(pt['holeVelocity_m_s']):6.2f} {float(pt['deltaP_Pa']):7.1f}  {pt['constraintState']}"
        )

    if optimization_list:
        lines.append("")
        lines.append("Optimization-ranked summer points:")
        lines.append("  Qtot(L/min)  Tmax(C)  Tmean(C)  dT(C)  cooling(W)  water(kg/d)  Pcool(W)  dP(Pa)  state")
        for pt in optimization_list:
            tmax = max(float(pt['TbottomFullPower_C']), float(pt['TtopFullPower_C']))
            tmean = 0.5 * (float(pt['TbottomFullPower_C']) + float(pt['TtopFullPower_C']))
            total_power = sum(
                value for value in (
                    float(pt.get("spotCoolerPower_W", np.nan)),
                    float(pt.get("assistBlowerPower_W", np.nan)),
                )
                if np.isfinite(value)
            )
            lines.append(
                f"  {float(pt['totalFlow_Lpm']):10.1f} {tmax:8.1f} {tmean:9.1f} "
                f"{abs(float(pt['TbottomFullPower_C']) - float(pt['TtopFullPower_C'])):6.1f} {max(-float(pt['QtoBed_W']), 0.0):10.1f} "
                f"{float(pt['waterLoss_kg_day']):11.2f} {total_power:9.1f} {float(pt['deltaP_Pa']):7.1f}  {pt['constraintState']}"
            )

    lines.append("============================================================")
    (output_dir / "study_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_summary(payload: dict, config: dict, output_dir: Path) -> None:
    rec = payload["recommended"]
    best_available = payload.get("best_available", {})
    reference = active_reference_point(payload)
    cfg = payload["active_configuration"]
    limits = effective_limits(payload, config)
    sweep = payload["sweep"]
    summary_inputs = config.get("summary_inputs", {})
    opt_cfg = optimization_config(config)
    model_overrides = normalized_model_overrides(config)
    env_inputs = effective_environment_inputs(summary_inputs, model_overrides)
    lumped_loss = dict(payload.get("lumped_loss", {}))
    bin_inputs = effective_bin_inputs(summary_inputs, model_overrides)
    aeration_inputs = effective_aeration_inputs(summary_inputs, model_overrides)
    heater_inputs = effective_heater_inputs(summary_inputs, model_overrides)
    wall_inputs = effective_bin_wall_inputs(summary_inputs, model_overrides)
    wire_inputs = summary_inputs.get("wire", {})
    effective_summary_inputs = dict(summary_inputs)
    effective_summary_inputs["bin"] = dict(bin_inputs)
    effective_summary_inputs["bin_wall"] = dict(wall_inputs)
    effective_summary_inputs["heater_tube"] = dict(heater_inputs)
    effective_summary_inputs["aeration"] = dict(aeration_inputs)
    perforation_inputs = model_overrides.get("perforation", {})
    evaporation_inputs = model_overrides.get("evaporation", {})
    model_aeration = model_overrides.get("aeration", {})
    wetted_area = derived_wetted_area_m2(summary_inputs, model_overrides)
    tube_layout = compute_tube_layout_and_habitat(payload, config)

    feasible_count = sum(1 for pt in rank_points_by_criteria(payload, config) if pt.get("selectionStage") == "feasible")
    top_n = int(opt_cfg.get("top_candidate_count", 8))
    top_list = matlab_top_candidates(payload, config, top_n)
    optimization_list = optimization_ranked_points(payload, config)
    consistency_warnings = export_consistency_warnings(payload, config)

    lines: list[str] = []
    lines.append("============================================================")
    lines.append("VERMICOMPOSTER AIR-HEATER DESIGN MODEL")
    lines.append("============================================================")
    append_line(lines, "Greenhouse air design temperature", f"{float(env_inputs['greenhouseAir_C']):.1f} C")
    append_line(lines, "Active vessel configuration", display_configuration_label(str(cfg["label"])))
    lines.append("")
    lines.append("Design specifications:")
    if {"length_m", "width_m", "height_m", "fillFraction"} <= bin_inputs.keys():
        lines.append(
            "  Vessel interior                         : "
            f"{float(bin_inputs['length_m']):.2f} m x {float(bin_inputs['width_m']):.2f} m x "
            f"{float(bin_inputs['height_m']):.2f} m, fill fraction {float(bin_inputs['fillFraction']):.2f}"
        )
    if not bool(bin_inputs.get("openTop", str(cfg["label"]).strip().lower() == "open top")):
        vent_count = int(bin_inputs.get("ventHoleCount", 0))
        vent_diameter_m = float(bin_inputs.get("ventHoleDiameter_m", 0.0))
        if vent_count > 0 and vent_diameter_m > 0:
            vent_area_per_hole = np.pi * vent_diameter_m**2 / 4.0
            lines.append(
                "  Vent specification                     : "
                f"{vent_count} vent hole(s), {1000.0 * vent_diameter_m:.1f} mm diameter each, "
                f"area per hole {vent_area_per_hole:.6f} m^2, total vent area {vent_count * vent_area_per_hole:.6f} m^2"
            )
    if {"sheetThickness_mm", "sheetK_W_mK", "insulationThickness_mm", "insulationK_W_mK"} <= wall_inputs.keys():
        lines.append(
            "  Vessel wall stack                       : "
            f"{float(wall_inputs['sheetThickness_mm']):.1f} mm steel (k = {float(wall_inputs['sheetK_W_mK']):.1f} W/m-K) + "
            f"{float(wall_inputs['insulationThickness_mm']):.1f} mm insulation (k = {float(wall_inputs['insulationK_W_mK']):.3f} W/m-K)"
        )
    if {"length_m", "ID_mm", "OD_mm", "insulationThickness_mm"} <= heater_inputs.keys():
        lines.append(
            "  Heater tube                             : "
            f"L {float(heater_inputs['length_m']):.3f} m, ID {float(heater_inputs['ID_mm']):.1f} mm, "
            f"OD {float(heater_inputs['OD_mm']):.1f} mm, insulation {float(heater_inputs['insulationThickness_mm']):.1f} mm"
        )
    if {"coilMeanD_mm", "coilSpan_m"} <= heater_inputs.keys():
        lines.append(
            "  Heater coil geometry                    : "
            f"mean D {float(heater_inputs['coilMeanD_mm']):.1f} mm, pitch {float(reference.get('pitch_mm', np.nan)):.1f} mm, "
            f"axial span {1000 * float(heater_inputs['coilSpan_m']):.1f} mm"
        )
    if {"nParallelTubes", "length_m", "ID_mm", "OD_mm"} <= aeration_inputs.keys():
        lines.append(
            "  Aeration network                        : "
            f"{int(aeration_inputs['nParallelTubes'])} tubes, each L {float(aeration_inputs['length_m']):.2f} m, "
            f"ID {float(aeration_inputs['ID_mm']):.1f} mm, OD {float(aeration_inputs['OD_mm']):.1f} mm"
        )
        header_enabled = bool(aeration_inputs.get("headerEnabled", True))
        splitter_inlet_mm = float(aeration_inputs.get("splitterInletID_mm", aeration_inputs.get("branchConnectorID_mm", 19.0)))
        connector_id_mm = float(aeration_inputs.get("branchConnectorID_mm", 19.0))
        connector_length_m = float(aeration_inputs.get("branchConnectorLength_m", 0.0))
        splitter_outlets = int(round(float(aeration_inputs.get("splitterOutletCount", 4))))
        if header_enabled:
            lines.append(
                "  Distribution manifold                  : "
                f"header ID {float(aeration_inputs.get('headerID_mm', np.nan)):.1f} mm feeding "
                f"{splitter_outlets}-way splitter(s), splitter inlet ID {splitter_inlet_mm:.1f} mm, "
                f"branch connector ID {connector_id_mm:.1f} mm, connector length {connector_length_m:.3f} m"
            )
        else:
            lines.append(
                "  Distribution manifold                  : "
                f"no separate header; {splitter_outlets}-way splitter(s), splitter inlet ID {splitter_inlet_mm:.1f} mm, "
                f"branch connector ID {connector_id_mm:.1f} mm, connector length {connector_length_m:.3f} m"
            )
        lines.append(
            "  Contraction model                      : "
            f"Idel'chik conical converging bellmouth without end wall, "
            f"{float(aeration_inputs.get('contractionConeAngle_deg', 60.0)):.0f} deg"
        )
    if tube_layout:
        center_pitch = tube_layout["pitches_m"]
        clear_gap = tube_layout["clear_gaps_m"]
        if center_pitch.size > 0:
            pitch_text = ", ".join(f"{1000.0 * value:.1f}" for value in center_pitch)
            gap_text = ", ".join(f"{1000.0 * value:.1f}" for value in clear_gap)
        else:
            pitch_text = "n/a"
            gap_text = "n/a"
        lines.append(
            "  Tube layout                             : "
            f"center height {1000.0 * float(tube_layout['y_centers_m'][0]):.1f} mm, "
            f"horizontal pitch {pitch_text} mm, clear OD gap {gap_text} mm, "
            f"nearest center spacing {1000.0 * float(tube_layout.get('nearest_center_spacing_m', np.nan)):.1f} mm"
        )
    if {"holesPerTube", "holeDiameter_m"} <= perforation_inputs.keys():
        hole_area_m2 = np.pi * float(perforation_inputs["holeDiameter_m"]) ** 2 / 4.0
        total_holes = int(aeration_inputs.get("nParallelTubes", 1)) * int(perforation_inputs["holesPerTube"])
        lines.append(
            "  Perforations                            : "
            f"{int(perforation_inputs['holesPerTube'])} holes/tube, {1000 * float(perforation_inputs['holeDiameter_m']):.1f} mm dia, "
            f"total holes {total_holes}, total open area {total_holes * hole_area_m2:.6f} m^2"
        )
    if evaporation_inputs:
        lines.append(
            "  Moisture model                          : "
            f"RH {100 * float(evaporation_inputs.get('relativeHumidity', 0.0)):.0f}%, "
            f"wetted area {wetted_area:.3f} m^2, "
            f"surface temp {float(evaporation_inputs.get('surfaceTemp_C', np.nan)):.1f} C"
        )
    if "bedH_W_m2K" in model_aeration:
        lines.append(
            "  Bed-side heat-transfer assumption       : "
            f"effective sensible h = {float(model_aeration['bedH_W_m2K']):.1f} W/m^2-K"
        )
    if "releaseFraction" in model_aeration:
        lines.append(
            "  Perforation release fraction            : "
            f"{100 * float(model_aeration['releaseFraction']):.1f}% of tube flow"
        )
    if {"name", "gauge", "diameter_mm"} <= wire_inputs.keys():
        lines.append(
            "  Heater wire                             : "
            f"{wire_inputs['name']}, {wire_inputs['gauge']}, {float(wire_inputs['diameter_mm']):.3f} mm diameter"
        )
    lines.append(
        "  Operating sweep                         : "
        f"V = {format_sweep(sweep['voltage_V'])} V, flow = {format_sweep(sweep['totalFlow_Lpm'])} L/min"
    )
    lines.append(
        "  Limits                                  : "
        f"I_total <= {float(limits['maxTotalCurrent_A']):.1f} A, "
        f"T_wire <= {float(limits['maxWireTemp_C']):.1f} C, "
        f"T_air,out <= {float(limits['maxAirOutletTemp_C']):.1f} C, "
        f"dP <= {float(limits['maxPressureDrop_Pa']):.0f} Pa, "
        f"water loss <= {float(limits.get('maxWaterLoss_kg_day', np.nan)):.1f} kg/24 h, "
        f"hole velocity <= {float(limits.get('maxHoleVelocity_m_s', np.inf)):.2f} m/s"
    )
    if "bedSetpoint_C" in summary_inputs:
        append_line(lines, f"Lumped heat-loss requirement at {float(summary_inputs['bedSetpoint_C']):.1f} C bed", f"{float(payload['requiredHeat_W']):.1f} W")
    else:
        append_line(lines, "Lumped heat-loss requirement", f"{float(payload['requiredHeat_W']):.1f} W")
    if "UA_W_K" in lumped_loss:
        append_line(lines, "Overall UA to greenhouse air", f"{float(lumped_loss['UA_W_K']):.3f} W/K")
    elif "overallUA_W_K" in summary_inputs:
        append_line(lines, "Overall UA to greenhouse air", f"{float(summary_inputs['overallUA_W_K']):.3f} W/K")
    if "externalH_W_m2K" in lumped_loss:
        append_line(lines, "External natural-convection h estimate", f"{float(lumped_loss['externalH_W_m2K']):.2f} W/m^2-K")
    elif "externalNaturalConvection_h_W_m2K" in summary_inputs:
        append_line(lines, "External natural-convection h estimate", f"{float(summary_inputs['externalNaturalConvection_h_W_m2K']):.2f} W/m^2-K")
    if "tau_h" in lumped_loss:
        append_line(lines, "Lumped time constant", f"{float(lumped_loss['tau_h']):.1f} h")
    elif "lumpedTimeConstant_h" in summary_inputs:
        append_line(lines, "Lumped time constant", f"{float(summary_inputs['lumpedTimeConstant_h']):.1f} h")
    if "Bi_lumped" in lumped_loss:
        append_line(lines, "Lumped Biot number estimate", f"{float(lumped_loss['Bi_lumped']):.2f}")
    elif "lumpedBiotNumber" in summary_inputs:
        append_line(lines, "Lumped Biot number estimate", f"{float(summary_inputs['lumpedBiotNumber']):.2f}")
    if "lumpedStrictlyValid" in lumped_loss:
        if bool(lumped_loss["lumpedStrictlyValid"]):
            append_line(lines, "Lumped-capacitance criterion", "Bi < 0.1 satisfied")
        else:
            append_line(lines, "Lumped-capacitance criterion", "Bi < 0.1 NOT satisfied")
            lines.append("                                             use the lumped model for")
            lines.append("                                             heat-loss sizing, not for")
            lines.append("                                             detailed internal gradients")
    elif "lumpedStrictlyValid" in summary_inputs:
        if bool(summary_inputs["lumpedStrictlyValid"]):
            append_line(lines, "Lumped-capacitance criterion", "Bi < 0.1 satisfied")
        else:
            append_line(lines, "Lumped-capacitance criterion", "Bi < 0.1 NOT satisfied")
            lines.append("                                             use the lumped model for")
            lines.append("                                             heat-loss sizing, not for")
            lines.append("                                             detailed internal gradients")
    append_line(
        lines,
        f"Minimum total flow at 100% duty, {float(limits['maxAirOutletTemp_C']):.1f} C cap",
        f"{float(payload['minimumFlowAt100Duty_Lpm']):.1f} L/min",
    )
    append_line(
        lines,
        f"Minimum total flow at {100 * float(sweep['targetDuty']):.0f}% duty, {float(limits['maxAirOutletTemp_C']):.1f} C cap",
        f"{float(payload['minimumFlowAtTargetDuty_Lpm']):.1f} L/min",
    )
    append_line(lines, "Maximum allowed total current", f"{float(limits['maxTotalCurrent_A']):.1f} A")
    append_line(lines, "Feasible operating points found", f"{feasible_count} of {len(payload['design_points'])}")
    if consistency_warnings:
        lines.append("")
        lines.append("Export/config consistency warnings:")
        for item in consistency_warnings:
            lines.append(f"  - {item}")

    lines.append("")
    lines.append("Vessel configuration comparison:")
    lines.append("  config                    UA(W/K)  Qreq(W)  Flow100  Flow20  Pref(W)  Qref(L/min)  Iref(A)  status")
    shown_points: list[tuple[float, float, float]] = []
    for item in payload["configuration_comparisons"]:
        rec_cfg = item.get("recommended", {})
        ref_cfg = item.get("best_available", rec_cfg)
        shown = rec_cfg if rec_cfg else ref_cfg
        shown_points.append(
            (
                float(shown["pitch_mm"]),
                float(shown["voltage_V"]),
                float(shown["totalFlow_Lpm"]),
            )
        )
        lines.append(
            f"  {item['label']:<24} {float(item['UA_W_K']):7.3f}  {float(item['Qreq_W']):7.1f}  "
            f"{float(item['flow100_Lpm']):7.1f} {float(item['flowTargetDuty_Lpm']):7.1f} "
            f"{float(shown['totalPower_W']):8.1f} {float(shown['totalFlow_Lpm']):11.1f} "
            f"{float(shown['totalCurrent_A']):8.2f}  {recommended_status_short(rec_cfg if rec_cfg else shown, limits, opt_cfg)}"
        )
    if shown_points and all(
        np.isclose(triple[0], shown_points[0][0])
        and np.isclose(triple[1], shown_points[0][1])
        and np.isclose(triple[2], shown_points[0][2])
        for triple in shown_points[1:]
    ):
        lines.append("  Note                               : all three comparison rows collapse to the same")
        lines.append("                                       best-available heater operating point because no")
        lines.append("                                       hard-safe configuration point exists in this sweep;")
        lines.append("                                       the fallback ranking is then dominated by the same")
        lines.append("                                       heater-side outlet/current/water-loss boundary")

    append_selection_process(lines, payload, config)
    append_optimization_process(lines, payload, config)
    append_moisture_process(lines, payload, config)
    append_radial_profile_process(lines, payload, config)
    append_tube_layout_process(lines, payload, config, tube_layout)
    append_geometry_trade_summary(lines, payload, config)

    write_design_point_block(
        lines,
        "Criteria-based recommended design point:",
        str(cfg["label"]),
        rec,
        effective_summary_inputs,
        limits,
        opt_cfg,
        sweep,
        configuration_extras(config, str(cfg["label"])),
    )
    write_design_point_block(
        lines,
        "Best-available active reference point:",
        str(cfg["label"]),
        best_available,
        effective_summary_inputs,
        limits,
        opt_cfg,
        sweep,
        configuration_extras(config, str(cfg["label"])),
    )

    vent_entry = next((item for item in payload["configuration_comparisons"] if item["label"].strip().lower() == "covered top + vent"), None)
    if vent_entry is not None:
        write_design_point_block(
            lines,
            "Covered-top-plus-vent criteria-based design point:",
            vent_entry["label"],
            vent_entry.get("recommended", {}),
            effective_summary_inputs,
            limits,
            opt_cfg,
            sweep,
            configuration_extras(config, vent_entry["label"]),
        )
        write_design_point_block(
            lines,
            "Covered-top-plus-vent best-available reference point:",
            vent_entry["label"],
            vent_entry.get("best_available", {}),
            effective_summary_inputs,
            limits,
            opt_cfg,
            sweep,
            configuration_extras(config, vent_entry["label"]),
        )

    lines.append("")
    if rec:
        lines.append("Top candidate points for active configuration:")
    else:
        lines.append("Top best-available candidate points for active configuration:")
    lines.append("  pitch  V   Qtot(L/min)  Ptot(W)  Vhole  evap(kg/d)  Tair(C)  Twire(C)  Qbed(W)  duty   dP(Pa)  stage")
    for pt in top_list:
        lines.append(
            f"  {float(pt['pitch_mm']):5.1f} {float(pt['voltage_V']):2.0f}    "
            f"{float(pt['totalFlow_Lpm']):8.1f}   {float(pt['totalPower_W']):7.1f}   {float(pt['holeVelocity_m_s']):5.2f}     {float(pt['waterLoss_kg_day']):6.2f}   "
            f"{float(pt['airOutlet_C']):7.1f}   {float(pt['wireMax_C']):8.1f}   {float(pt['QtoBed_W']):7.1f}  "
            f"{float(pt['dutyNeeded']):4.2f}   {float(pt['deltaP_Pa']):6.1f}  {pt['selectionStage']}"
        )

    if optimization_list:
        lines.append("")
        lines.append("Optimization-ranked points:")
        lines.append("  pitch  V   Qtot(L/min)  Ptot(W)  Tb(C)  Tt(C)  dT(C)  Vhole  evap(kg/d)  Tair(C)  Twire(C)  duty  dP(Pa)  state")
        for pt in optimization_list:
            lines.append(
                f"  {float(pt['pitch_mm']):5.1f} {float(pt['voltage_V']):3.0f}  "
                f"{float(pt['totalFlow_Lpm']):10.1f} {float(pt['totalPower_W']):8.1f} "
                f"{float(pt.get('TbottomFullPower_C', np.nan)):6.1f} {float(pt.get('TtopFullPower_C', np.nan)):6.1f} "
                f"{float(pt['tempSpread_C']):6.1f} {float(pt['holeVelocity_m_s']):5.2f}   {float(pt['waterLoss_kg_day']):8.2f} {float(pt['airOutlet_C']):8.1f} "
                f"{float(pt['wireMax_C']):9.1f} {float(pt['dutyNeeded']):5.2f} {float(pt['deltaP_Pa']):7.1f}  "
                f"{pt['constraintState']}"
            )

    lines.append("")
    lines.append("Constraint logic:")
    lines.append("  Hard-safe means all of the following are satisfied simultaneously:")
    lines.append(
        "  Explicit equipment limits            : "
        f"I_total <= {float(limits['maxTotalCurrent_A']):.1f} A, "
        f"T_air,out <= {float(limits['maxAirOutletTemp_C']):.1f} C, "
        f"T_wire <= {float(limits['maxWireTemp_C']):.1f} C, "
        f"dP <= {float(limits['maxPressureDrop_Pa']):.0f} Pa, "
        f"water loss <= {float(limits.get('maxWaterLoss_kg_day', np.nan)):.1f} kg/24 h, "
        f"hole velocity <= {float(limits.get('maxHoleVelocity_m_s', np.inf)):.2f} m/s"
    )
    lines.append(
        "  Explicit thermal limits              : "
        f"T_bottom >= {float(opt_cfg['min_bottom_temp_C']):.1f} C, "
        f"T_top >= {float(opt_cfg['min_top_temp_C']):.1f} C, "
        f"|T_bottom - T_top| <= {float(opt_cfg['spread_limit_C']):.1f} C, "
        "Q_to_bed > 0"
    )
    lines.append("  Feasible means hard-safe plus meeting the full heating target:")
    lines.append("                                         Q_to_bed >= lumped heat-loss requirement")
    lines.append("                                         equivalently, dutyNeeded <= 1.00 at full power")
    lines.append("  Recommended design point             : best lexicographic point among hard-safe candidates only")
    lines.append("  Best-available reference point       : best lexicographic point in the full sweep when no hard-safe")
    lines.append("                                         candidate exists; used for plots and fallback diagnostics")
    lines.append("  Clean ordering inside the hard-safe set")
    lines.append("                                         : maximize min(T_bottom, T_top), then minimize power,")
    lines.append("                                           then minimize |T_bottom - T_top|, then break remaining")
    lines.append("                                           ties with mean bed temperature and heat shortfall")
    lines.append("  Failure categorization order         : current, wire temperature, outlet air, pressure drop,")
    lines.append("                                         water loss, hole velocity, bottom temperature, top temperature,")
    lines.append("                                         spread, then residual heat shortfall/nonpositive heat")
    lines.append("============================================================")

    (output_dir / "study_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_root_summary_index(root_output_dir: Path, mode_outputs: dict[str, Path]) -> None:
    lines = [
        "============================================================",
        "HEAT TRANSFER STUDY OUTPUT INDEX",
        "============================================================",
        "This root summary points to the mode-specific output folders.",
        "",
    ]

    if "heating" in mode_outputs:
        lines.append("Heating mode:")
        lines.append(f"  Folder                             : {mode_outputs['heating']}")
        lines.append(f"  Summary                            : {mode_outputs['heating'] / 'study_summary.txt'}")
        lines.append("")
    if "cooling" in mode_outputs:
        lines.append("Cooling mode:")
        lines.append(f"  Folder                             : {mode_outputs['cooling']}")
        lines.append(f"  Summary                            : {mode_outputs['cooling'] / 'study_summary.txt'}")
        lines.append("")
    if "year_round" in mode_outputs:
        lines.append("Year-round analysis:")
        lines.append(f"  Folder                             : {mode_outputs['year_round']}")
        lines.append(f"  Summary                            : {mode_outputs['year_round'] / 'study_summary.txt'}")
        lines.append("")

    if not mode_outputs:
        lines.append("No run modes were enabled in study_config.json.")

    (root_output_dir / "study_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    config_path = args.config.resolve()
    config = load_json(config_path)
    root_output_dir = resolve_path(config_path.parent, config.get("output_dir", "outputs"))
    root_output_dir.mkdir(parents=True, exist_ok=True)
    run_modes = enabled_run_modes(config)
    mode_outputs: dict[str, Path] = {}
    heating_payload: dict | None = None
    cooling_payload: dict | None = None

    if run_modes.get("heating", False):
        output_dir = root_output_dir / "heating"
        output_dir.mkdir(parents=True, exist_ok=True)
        auto_refresh_export_if_needed(config_path, config)
        data_json = resolve_path(config_path.parent, config["data_json"])
        grid_points = int(config.get("plot_points", 300))
        curve_count = int(config.get("constraint_curve_count", 10))
        wire_plot_points = int(config.get("wire_plot_points", 100))

        payload = load_json(data_json)
        payload = apply_config_overrides(payload, config)
        active_selected = select_recommended_point(payload, config)
        payload["recommended"] = active_selected
        payload["best_available"] = select_best_available_point(payload, config)
        for item in payload.get("configuration_comparisons", []):
            pts = item.get("design_points", [])
            selected = select_recommended_point(payload, config, pts)
            item["recommended"] = selected
            item["best_available"] = select_best_available_point(payload, config, pts)

        (output_dir / "vessel_configuration_comparison.png").unlink(missing_ok=True)
        (output_dir / "summer_spot_cooler_curves.png").unlink(missing_ok=True)
        (output_dir / "summer_evaporative_cooling_curves.png").unlink(missing_ok=True)
        plot_constraint_maps(payload, output_dir, grid_points, curve_count)
        plot_total_current_curves(payload, output_dir, curve_count)
        plot_moisture_maps(payload, output_dir, grid_points, curve_count)
        plot_radial_air_cooling_profile(payload, config, output_dir)
        plot_habitat_exclusion_cross_section(payload, config, output_dir)
        plot_wire_trade_study(
            payload,
            config,
            output_dir,
            wire_plot_points,
            "wire_diameter_sweep",
            "Wire Diameter Trade Study",
            "wire_diameter_trade_study.png",
        )
        plot_wire_trade_study(
            payload,
            config,
            output_dir,
            wire_plot_points,
            "coil_diameter_sweep",
            "Coil Mean Diameter Trade Study",
            "coil_diameter_trade_study.png",
        )
        plot_wire_trade_study(
            payload,
            config,
            output_dir,
            wire_plot_points,
            "wire_length_sweep",
            "Wire Length Trade Study",
            "wire_length_trade_study.png",
        )
        write_summary(payload, config, output_dir)
        mode_outputs["heating"] = output_dir
        heating_payload = payload

    if run_modes.get("cooling", False):
        output_dir = root_output_dir / "cooling"
        output_dir.mkdir(parents=True, exist_ok=True)
        payload = build_summer_cooling_payload(config)
        for stale in (
            "constraint_performance_maps.png",
            "total_current_curves.png",
            "moisture_perforation_maps.png",
            "wire_diameter_trade_study.png",
            "coil_diameter_trade_study.png",
            "wire_length_trade_study.png",
            "vessel_configuration_comparison.png",
            "summer_spot_cooler_curves.png",
            "summer_evaporative_cooling_curves.png",
        ):
            (output_dir / stale).unlink(missing_ok=True)
        plot_summer_cooling_curves(payload, config, output_dir)
        plot_radial_air_cooling_profile(payload, config, output_dir)
        plot_habitat_exclusion_cross_section(payload, config, output_dir)
        write_cooling_summary(payload, config, output_dir)
        mode_outputs["cooling"] = output_dir
        cooling_payload = payload

    if year_round_analysis_config(config).get("enabled", False) and heating_payload is not None and cooling_payload is not None:
        output_dir = root_output_dir / "year_round"
        output_dir.mkdir(parents=True, exist_ok=True)
        annual_payload = build_year_round_payload(heating_payload, cooling_payload, config)
        plot_year_round_energy_cost(annual_payload, output_dir)
        write_year_round_summary(annual_payload, output_dir)
        mode_outputs["year_round"] = output_dir

    write_root_summary_index(root_output_dir, mode_outputs)
    print(f"Wrote standalone heat-transfer study outputs to {root_output_dir}")


if __name__ == "__main__":
    main()

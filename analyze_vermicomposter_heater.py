from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d


PROJECT_DIR = Path(__file__).resolve().parent
MODEL_NAME = "vermicomposter_air_heater_design_model"
DEFAULT_OUTPUT_DIR = PROJECT_DIR / "python_heater_analysis_moisture"
DEFAULT_STUDY_CONFIG = PROJECT_DIR / "Heat Transfer Study" / "study_config.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Export results from the MATLAB vermicomposter heater model and "
            "regenerate the key plots in Python with clearer recommendation "
            "definitions and 300-point display resolution."
        )
    )
    parser.add_argument(
        "--plot-points",
        type=int,
        default=300,
        help=(
            "Display resolution used for continuous Python plots. "
            "This is interpolation density, not a new physics solve."
        ),
    )
    parser.add_argument(
        "--model-grid-points",
        type=int,
        default=None,
        help=(
            "Optional override for MATLAB sweep density on voltage, flow, and "
            "wire diameter. This can be very slow for large values."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=(
            "Retained for backward compatibility, but outputs are always written "
            "to the fixed tracked analysis directory."
        ),
    )
    parser.add_argument(
        "--matlab",
        type=Path,
        default=None,
        help="Optional explicit path to matlab.exe.",
    )
    parser.add_argument(
        "--study-config",
        type=Path,
        default=DEFAULT_STUDY_CONFIG,
        help="Optional shared JSON config file used to populate MATLAB overrides.",
    )
    return parser.parse_args()


def find_matlab(explicit_path: Path | None) -> Path:
    if explicit_path is not None:
        return explicit_path

    path_match = shutil.which("matlab")
    if path_match:
        return Path(path_match)

    fallbacks = [
        Path(r"C:\Program Files\MATLAB\R2025b\bin\matlab.exe"),
        Path(r"C:\Program Files\MATLAB\R2025a\bin\matlab.exe"),
    ]
    for candidate in fallbacks:
        if candidate.exists():
            return candidate

    raise FileNotFoundError("Unable to locate matlab.exe.")


def matlab_path_string(path: Path) -> str:
    return path.as_posix().replace("'", "''")


def load_json(path: Path | None) -> dict:
    if path is None or not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def resolve_path(root: Path, raw_path: str) -> Path:
    path = Path(raw_path)
    if path.is_absolute():
        return path
    return (root / path).resolve()


def reset_generated_output_dir(output_dir: Path) -> None:
    preserved_names = {".gitattributes", ".gitignore", ".gitkeep"}
    output_dir.mkdir(parents=True, exist_ok=True)

    for child in output_dir.iterdir():
        if child.name in preserved_names:
            continue
        if child.is_dir():
            shutil.rmtree(child)
        else:
            child.unlink()


def model_overrides_signature(overrides: dict) -> str:
    encoded = json.dumps(overrides, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


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


def normalize_model_overrides_for_matlab(overrides: dict) -> dict:
    normalized = apply_override_aliases(strip_comment_keys(json.loads(json.dumps(overrides))))
    sweep_cfg = normalized.get("sweep")
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
    return normalized


def normalize_matlab_parallel_config(parallel_cfg: dict) -> dict:
    normalized = strip_comment_keys(json.loads(json.dumps(parallel_cfg or {})))
    if not isinstance(normalized, dict):
        return {}
    out = dict(normalized)
    if "enabled" in out:
        out["enabled"] = bool(out["enabled"])
    if "workers" in out and out["workers"] not in (None, ""):
        out["workers"] = int(out["workers"])
    if "autoStartPool" in out:
        out["autoStartPool"] = bool(out["autoStartPool"])
    loops = dict(out.get("loops", {}))
    loops.setdefault("exportSweepRowsParfor", bool(out.get("enabled", False)))
    out["loops"] = loops
    out["enabled"] = bool(out.get("enabled", False)) and bool(loops.get("exportSweepRowsParfor", True))
    return out


def matlab_parallel_config(study_config: dict) -> dict:
    grouped = study_config.get("parallel", {})
    if isinstance(grouped, dict):
        raw = grouped.get("matlab", study_config.get("matlab_parallel", {}))
    else:
        raw = study_config.get("matlab_parallel", {})
    return normalize_matlab_parallel_config(raw)


def matlab_literal(value) -> str:
    if isinstance(value, dict):
        if not value:
            return "struct()"
        parts = []
        for key, item in value.items():
            parts.append(f"'{key}', {matlab_literal(item)}")
        return "struct(" + ", ".join(parts) + ")"
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return "[]"
    if isinstance(value, str):
        return "'" + value.replace("'", "''") + "'"
    if isinstance(value, (int, float)):
        return repr(value)
    if isinstance(value, list):
        if not value:
            return "[]"
        if all(isinstance(item, (int, float)) for item in value):
            return "[" + " ".join(repr(item) for item in value) + "]"
        return "{" + ", ".join(matlab_literal(item) for item in value) + "}"
    raise TypeError(f"Unsupported override type: {type(value)!r}")


def build_export_script(
    export_json: Path,
    model_grid_points: int | None,
    helper_name: str,
    study_overrides: dict,
    matlab_parallel: dict,
) -> str:
    export_json_str = matlab_path_string(export_json)
    project_dir_str = matlab_path_string(PROJECT_DIR)
    sweep_cfg = study_overrides.get("sweep", {}) if isinstance(study_overrides, dict) else {}
    voltage_values = sweep_cfg.get("voltage_V") if isinstance(sweep_cfg, dict) else None
    flow_values = sweep_cfg.get("totalFlow_Lpm") if isinstance(sweep_cfg, dict) else None
    wire_values = sweep_cfg.get("wireDiameter_mm") if isinstance(sweep_cfg, dict) else None
    v_min = float(voltage_values[0]) if isinstance(voltage_values, list) and voltage_values else 2.0
    v_max = float(voltage_values[-1]) if isinstance(voltage_values, list) and voltage_values else 120.0
    q_min = float(flow_values[0]) if isinstance(flow_values, list) and flow_values else 40.0
    q_max = float(flow_values[-1]) if isinstance(flow_values, list) and flow_values else 400.0
    d_min = float(wire_values[0]) if isinstance(wire_values, list) and wire_values else 0.18
    d_max = float(wire_values[-1]) if isinstance(wire_values, list) and wire_values else 0.64
    if model_grid_points is None:
        grid_override = "modelGridPoints = [];"
    else:
        grid_override = f"modelGridPoints = {int(model_grid_points)};"
    override_literal = matlab_literal(study_overrides) if study_overrides else "struct()"
    parallel_literal = matlab_literal(matlab_parallel) if matlab_parallel else "struct()"

    return textwrap.dedent(
        f"""
        function {helper_name}
        addpath('{project_dir_str}');
        exportJsonPath = '{export_json_str}';
        {grid_override}

        overrides = {override_literal};
        parallelCfg = {parallel_literal};
        if ~isfield(overrides, 'options') || ~isstruct(overrides.options)
            overrides.options = struct();
        end
        overrides.options.makePlots = false;
        overrides.options.includeTransient = false;
        overrides.options.printSummary = false;
        if ~isempty(fieldnames(parallelCfg))
            overrides.parallel = parallelCfg;
        end

        if ~isempty(modelGridPoints)
            if ~isfield(overrides, 'sweep') || ~isstruct(overrides.sweep)
                overrides.sweep = struct();
            end
            overrides.sweep.voltage_V = linspace({v_min}, {v_max}, modelGridPoints);
            overrides.sweep.totalFlow_Lpm = linspace({q_min}, {q_max}, modelGridPoints);
            overrides.sweep.wireDiameter_mm = linspace({d_min}, {d_max}, modelGridPoints);
        end

        results = {MODEL_NAME}(overrides);

        activeIdx = results.activeConfigurationIndex;
        if isfinite(activeIdx) && activeIdx >= 1 && ...
                activeIdx <= numel(results.configurationComparisons)
            activeCfg = results.configurationComparisons(activeIdx);
            activeScenario = activeCfg.scenario;
            rec = activeScenario.recommended;
        else
            activeCfg = struct();
            activeCfg.label = results.binModel.topMode;
            activeCfg.color = [0.15 0.15 0.15];
            activeScenario = results;
            rec = results.recommended;
            activeIdx = NaN;
        end

        payload = struct();
        payload.meta = struct();
        payload.meta.model_name = '{MODEL_NAME}';
        payload.meta.model_grid_points = modelGridPoints;
        payload.meta.plot_note = ['Python figures use interpolated display grids; ' ...
            'they do not change the underlying MATLAB physics unless model_grid_points is set.'];
        payload.meta.recommended_definition = ['Recommended = the single hard-constraint-safe operating point ' ...
            'selected for a vessel configuration from the MATLAB sweep by the shared lexicographic criteria ' ...
            'used throughout the study. Inside the hard-safe set, the rule is: maximize the colder-node ' ...
            'bed temperature, then minimize power, then minimize spread. If no hard-safe point exists, ' ...
            'the best-available point is retained only as a reference point.'];
        payload.meta.star_definition = ['The standalone Heat Transfer Study plots do not place stars ' ...
            'on the contour panels. Recommended operating points are reported in the text summary.'];

        payload.active_configuration = struct();
        payload.active_configuration.index = activeIdx;
        payload.active_configuration.label = activeCfg.label;
        payload.active_configuration.color = activeCfg.color;

        payload.recommended = simplifyRecommended(activeScenario.recommended);
        payload.best_available = simplifyRecommended(activeScenario.bestAvailable);
        payload.limits = activeScenario.params.limits;
        payload.sweep = activeScenario.params.sweep;
        payload.requiredHeat_W = activeScenario.requiredHeat_W;
        payload.minimumFlowAt100Duty_Lpm = activeScenario.minimumFlowAt100Duty_Lpm;
        payload.minimumFlowAtTargetDuty_Lpm = activeScenario.minimumFlowAtTargetDuty_Lpm;
        payload.lumped_loss = simplifyLumpedLoss(activeScenario.lumpedLoss, activeScenario.binModel);

        payload.configuration_comparisons = repmat(struct(), numel(results.configurationComparisons), 1);
        for i = 1:numel(results.configurationComparisons)
            cfg = results.configurationComparisons(i);
            payload.configuration_comparisons(i).label = cfg.label;
            payload.configuration_comparisons(i).color = cfg.color;
            payload.configuration_comparisons(i).UA_W_K = cfg.scenario.lumpedLoss.UA_W_K;
            payload.configuration_comparisons(i).Qreq_W = cfg.scenario.requiredHeat_W;
            payload.configuration_comparisons(i).flow100_Lpm = cfg.scenario.minimumFlowAt100Duty_Lpm;
            payload.configuration_comparisons(i).flowTargetDuty_Lpm = cfg.scenario.minimumFlowAtTargetDuty_Lpm;
            payload.configuration_comparisons(i).recommended = simplifyRecommended(cfg.scenario.recommended);
            payload.configuration_comparisons(i).best_available = simplifyRecommended(cfg.scenario.bestAvailable);
            payload.configuration_comparisons(i).status = shortRecommendedStatus( ...
                cfg.scenario.recommended, cfg.scenario.bestAvailable);
            payload.configuration_comparisons(i).design_points = simplifyDesignPoints(cfg.scenario.designPoints);
        end

        payload.design_points = simplifyDesignPoints(activeScenario.designPoints);
        payload.wire_diameter_sweep = simplifyTradeSweep(results.wireDiameterSweep, 'diameter_mm');
        payload.coil_diameter_sweep = simplifyTradeSweep(results.coilDiameterSweep, 'coilMeanD_mm');
        payload.wire_length_sweep = simplifyTradeSweep(results.wireLengthSweep, 'wireLength_m');

        jsonText = jsonencode(payload);
        fid = fopen(exportJsonPath, 'w');
        fwrite(fid, jsonText, 'char');
        fclose(fid);
        end

        function out = simplifyRecommended(rec)
            out = struct();
            if isempty(fieldnames(rec))
                return;
            end
            out.pitch_mm = rec.pitch_mm;
            out.voltage_V = rec.voltage_V;
            out.totalFlow_Lpm = rec.totalFlow_Lpm;
            out.totalPower_W = rec.totalPower_W;
            out.totalCurrent_A = rec.totalCurrent_A;
            out.wireLength_m = rec.wireLength_m;
            out.hWire_W_m2K = rec.hWire_W_m2K;
            out.airOutlet_C = rec.airOutlet_C;
            out.wireMax_C = rec.wireMax_C;
            out.deltaP_Pa = rec.deltaP_Pa;
            if isfield(rec, 'distributionDeltaP_Pa')
                out.distributionDeltaP_Pa = rec.distributionDeltaP_Pa;
                out.distributionHeader_Pa = rec.distributionHeader_Pa;
                out.distributionHeaderToSplitter_Pa = rec.distributionHeaderToSplitter_Pa;
                out.distributionSplitterBody_Pa = rec.distributionSplitterBody_Pa;
                out.distributionConnectorFriction_Pa = rec.distributionConnectorFriction_Pa;
                out.distributionConnectorToBranch_Pa = rec.distributionConnectorToBranch_Pa;
            end
            out.QtoBed_W = rec.QtoBed_W;
            out.dutyNeeded = rec.dutyNeeded;
            if isfield(rec, 'flowPerHole_Lpm')
                out.flowPerHole_Lpm = rec.flowPerHole_Lpm;
            end
            if isfield(rec, 'holeVelocity_m_s')
                out.holeVelocity_m_s = rec.holeVelocity_m_s;
            end
            if isfield(rec, 'waterLoss_kg_day')
                out.waterLoss_kg_day = rec.waterLoss_kg_day;
            end
            if isfield(rec, 'latentEvap_W')
                out.latentEvap_W = rec.latentEvap_W;
            end
            if isfield(rec, 'TbottomFullPower_C')
                out.TbottomFullPower_C = rec.TbottomFullPower_C;
            end
            if isfield(rec, 'TtopFullPower_C')
                out.TtopFullPower_C = rec.TtopFullPower_C;
            end
            out.isFeasible = rec.isFeasible;
            if isfield(rec, 'isConstraintSafe')
                out.isConstraintSafe = rec.isConstraintSafe;
            else
                out.isConstraintSafe = false;
            end
            if isfield(rec, 'meetsCurrentCap')
                out.meetsCurrentCap = rec.meetsCurrentCap;
            else
                out.meetsCurrentCap = false;
            end
        end

        function out = simplifyLumpedLoss(loss, binModel)
            out = struct();
            if isempty(fieldnames(loss))
                return;
            end
            out.requiredHeat_W = loss.requiredHeat_W;
            out.UA_W_K = loss.UA_W_K;
            if isfield(loss, 'tau_h')
                out.tau_h = loss.tau_h;
            end
            if isfield(loss, 'Bi_lumped')
                out.Bi_lumped = loss.Bi_lumped;
            end
            if isfield(loss, 'lumpedStrictlyValid')
                out.lumpedStrictlyValid = loss.lumpedStrictlyValid;
            end
            if nargin >= 2 && isstruct(binModel) && isfield(binModel, 'externalH_W_m2K')
                out.externalH_W_m2K = binModel.externalH_W_m2K;
            end
        end

        function out = simplifyDesignPoints(pts)
            out = repmat(struct(), numel(pts), 1);
            for j = 1:numel(pts)
                out(j).pitch_mm = pts(j).pitch_mm;
                out(j).nParallelTubes = pts(j).nParallelTubes;
                out(j).voltage_V = pts(j).voltage_V;
                out(j).totalFlow_Lpm = pts(j).totalFlow_Lpm;
                out(j).branchFlow_Lpm = pts(j).branchFlow_Lpm;
                out(j).totalPower_W = pts(j).totalPower_W;
                out(j).totalCurrent_A = pts(j).totalCurrent_A;
                out(j).branchResistance_Ohm = pts(j).branchResistance_Ohm;
                out(j).supplyEquivalentResistance_Ohm = pts(j).supplyEquivalentResistance_Ohm;
                out(j).branchCurrent_A = pts(j).branchCurrent_A;
                out(j).branchPower_W = pts(j).branchPower_W;
                if isfield(pts(j), 'parallelStringCount')
                    out(j).parallelStringCount = pts(j).parallelStringCount;
                end
                if isfield(pts(j), 'electricalTopologyLabel')
                    out(j).electricalTopologyLabel = pts(j).electricalTopologyLabel;
                end
                if isfield(pts(j), 'stringTubeCounts')
                    out(j).stringTubeCounts = pts(j).stringTubeCounts;
                end
                if isfield(pts(j), 'tubeCurrentMin_A')
                    out(j).tubeCurrentMin_A = pts(j).tubeCurrentMin_A;
                    out(j).tubeCurrentMean_A = pts(j).tubeCurrentMean_A;
                    out(j).tubeCurrentMax_A = pts(j).tubeCurrentMax_A;
                end
                if isfield(pts(j), 'tubePowerMin_W')
                    out(j).tubePowerMin_W = pts(j).tubePowerMin_W;
                    out(j).tubePowerMean_W = pts(j).tubePowerMean_W;
                    out(j).tubePowerMax_W = pts(j).tubePowerMax_W;
                end
                out(j).wireLength_m = pts(j).wireLength_m;
                out(j).hWire_W_m2K = pts(j).hWire_W_m2K;
                out(j).airOutlet_C = pts(j).airOutlet_C;
                out(j).wireMax_C = pts(j).wireMax_C;
                out(j).deltaP_Pa = pts(j).deltaP_Pa;
                if isfield(pts(j), 'distributionDeltaP_Pa')
                    out(j).distributionDeltaP_Pa = pts(j).distributionDeltaP_Pa;
                    out(j).distributionHeader_Pa = pts(j).distributionHeader_Pa;
                    out(j).distributionHeaderToSplitter_Pa = pts(j).distributionHeaderToSplitter_Pa;
                    out(j).distributionSplitterBody_Pa = pts(j).distributionSplitterBody_Pa;
                    out(j).distributionConnectorFriction_Pa = pts(j).distributionConnectorFriction_Pa;
                    out(j).distributionConnectorToBranch_Pa = pts(j).distributionConnectorToBranch_Pa;
                end
                out(j).branchQtoBed_W = pts(j).branchQtoBed_W;
                out(j).QtoBed_W = pts(j).QtoBed_W;
                out(j).dutyNeeded = pts(j).dutyNeeded;
                out(j).flowPerHole_Lpm = pts(j).flowPerHole_Lpm;
                out(j).holeVelocity_m_s = pts(j).holeVelocity_m_s;
                out(j).waterLoss_kg_day = pts(j).waterLoss_kg_day;
                out(j).latentEvap_W = pts(j).latentEvap_W;
                out(j).TbottomFullPower_C = pts(j).TbottomFullPower_C;
                out(j).TtopFullPower_C = pts(j).TtopFullPower_C;
                out(j).isFeasible = pts(j).isFeasible;
                out(j).isConstraintSafe = pts(j).isConstraintSafe;
            end
        end

        function out = simplifyTradeSweep(wireSweep, axisField)
            if isempty(wireSweep)
                out = struct([]);
                return;
            end
            out = repmat(struct(), numel(wireSweep), 1);
            for j = 1:numel(wireSweep)
                out(j).studyType = wireSweep(j).studyType;
                out(j).xField = wireSweep(j).xField;
                out(j).xLabel = wireSweep(j).xLabel;
                out(j).xUnit = wireSweep(j).xUnit;
                out(j).xValues = wireSweep(j).xValues;
                out(j).label = wireSweep(j).label;
                out(j).color = wireSweep(j).color;
                out(j).(axisField) = wireSweep(j).(axisField);
                out(j).targetPower_W = wireSweep(j).targetPower_W;
                out(j).targetFlow_Lpm = wireSweep(j).targetFlow_Lpm;
                out(j).targetCurrent_A = wireSweep(j).targetCurrent_A;
                out(j).optimalPower_W = wireSweep(j).optimalPower_W;
                out(j).optimalFlow_Lpm = wireSweep(j).optimalFlow_Lpm;
                out(j).optimalCurrent_A = wireSweep(j).optimalCurrent_A;
                out(j).optimalVoltage_V = wireSweep(j).optimalVoltage_V;
                out(j).optimalPitch_mm = wireSweep(j).optimalPitch_mm;
                out(j).optimalDuty = wireSweep(j).optimalDuty;
                out(j).optimalQtoBed_W = wireSweep(j).optimalQtoBed_W;
                out(j).optimalWireMax_C = wireSweep(j).optimalWireMax_C;
                out(j).optimalColderNode_C = wireSweep(j).optimalColderNode_C;
                out(j).selectedPower_W = wireSweep(j).selectedPower_W;
                out(j).selectedFlow_Lpm = wireSweep(j).selectedFlow_Lpm;
                out(j).selectedCurrent_A = wireSweep(j).selectedCurrent_A;
                out(j).selectedVoltage_V = wireSweep(j).selectedVoltage_V;
                out(j).selectedPitch_mm = wireSweep(j).selectedPitch_mm;
                out(j).selectedDuty = wireSweep(j).selectedDuty;
                out(j).selectedQtoBed_W = wireSweep(j).selectedQtoBed_W;
                out(j).selectedWireMax_C = wireSweep(j).selectedWireMax_C;
                out(j).selectedColderNode_C = wireSweep(j).selectedColderNode_C;
                out(j).selectedSelectionScore = wireSweep(j).selectedSelectionScore;
                out(j).selectionMode = wireSweep(j).selectionMode;
                out(j).bestAvailablePower_W = wireSweep(j).bestAvailablePower_W;
                out(j).bestAvailableFlow_Lpm = wireSweep(j).bestAvailableFlow_Lpm;
                out(j).bestAvailableCurrent_A = wireSweep(j).bestAvailableCurrent_A;
                out(j).bestAvailableVoltage_V = wireSweep(j).bestAvailableVoltage_V;
                out(j).bestAvailablePitch_mm = wireSweep(j).bestAvailablePitch_mm;
                out(j).bestAvailableDuty = wireSweep(j).bestAvailableDuty;
                out(j).bestAvailableQtoBed_W = wireSweep(j).bestAvailableQtoBed_W;
                out(j).bestAvailableWireMax_C = wireSweep(j).bestAvailableWireMax_C;
                out(j).bestAvailableColderNode_C = wireSweep(j).bestAvailableColderNode_C;
                out(j).bestAvailableMode = wireSweep(j).bestAvailableMode;
                out(j).bestIndex = wireSweep(j).bestIndex;
                out(j).bestAvailableIndex = wireSweep(j).bestAvailableIndex;
                out(j).isFeasible = wireSweep(j).isFeasible;
                out(j).isConstraintSafe = wireSweep(j).isConstraintSafe;
            end
        end

        function txt = shortRecommendedStatus(rec, bestAvailable)
            if ~isempty(fieldnames(rec))
                if rec.isFeasible
                    txt = 'feasible';
                elseif isfield(rec, 'isConstraintSafe') && rec.isConstraintSafe
                    txt = 'hard-safe';
                else
                    txt = 'best-available only';
                end
            elseif ~isempty(fieldnames(bestAvailable))
                txt = 'best-available only';
            else
                txt = 'none';
            end
        end
        """
    ).strip()


def export_results_json(
    matlab_exe: Path,
    output_dir: Path,
    model_grid_points: int | None,
    study_overrides: dict,
    matlab_parallel: dict,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    export_json = output_dir / "model_export.json"
    temp_export_dir = Path(tempfile.mkdtemp(prefix="heater_export_"))
    temp_export_json = temp_export_dir / "model_export.json"
    matlab_prefdir = temp_export_dir / "matlab_prefs"
    matlab_prefdir.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".m",
        prefix="pyexport_",
        delete=False,
    ) as tmp:
        helper_name = Path(tmp.name).stem
        script_text = build_export_script(
            temp_export_json,
            model_grid_points,
            helper_name,
            study_overrides,
            matlab_parallel,
        )
        tmp.write(script_text)
        temp_script = Path(tmp.name)

    batch_cmd = helper_name
    try:
        env = os.environ.copy()
        env["MATLAB_PREFDIR"] = str(matlab_prefdir)
        env["TMP"] = str(temp_export_dir)
        env["TEMP"] = str(temp_export_dir)
        subprocess.run(
            [str(matlab_exe), "-batch", batch_cmd],
            cwd=temp_script.parent,
            check=True,
            text=True,
            env=env,
        )
        if not temp_export_json.exists():
            raise FileNotFoundError(
                f"MATLAB export completed without creating the expected JSON: {temp_export_json}"
            )
        reset_generated_output_dir(output_dir)
        shutil.copy2(temp_export_json, export_json)
    finally:
        temp_script.unlink(missing_ok=True)
        shutil.rmtree(temp_export_dir, ignore_errors=True)

    return export_json


def load_payload(export_json: Path) -> dict:
    with export_json.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def write_payload(export_json: Path, payload: dict) -> None:
    export_json.write_text(json.dumps(payload), encoding="utf-8")


def annotate_export_payload(export_json: Path, normalized_overrides: dict, study_config_path: Path | None) -> None:
    payload = load_payload(export_json)
    meta = payload.setdefault("meta", {})
    meta["model_overrides_snapshot"] = normalized_overrides
    meta["model_overrides_signature"] = model_overrides_signature(normalized_overrides)
    if study_config_path is not None:
        meta["study_config_path"] = str(study_config_path)
    write_payload(export_json, payload)


def sync_export_to_study_data(export_json: Path, study_config_path: Path | None, study_config: dict) -> Path | None:
    if study_config_path is None or not study_config or "data_json" not in study_config:
        return None
    target = resolve_path(study_config_path.parent, str(study_config["data_json"]))
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(export_json, target)
    return target


def as_float_array(values) -> np.ndarray:
    return np.array([np.nan if value is None else float(value) for value in values], dtype=float)


def active_reference_point(payload: dict) -> dict:
    rec = payload.get("recommended", {})
    if isinstance(rec, dict) and rec.get("pitch_mm") is not None:
        return rec
    best = payload.get("best_available", {})
    if isinstance(best, dict) and best.get("pitch_mm") is not None:
        return best
    return {}


def resample_regular_grid(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, nx: int, ny: int, nearest: bool = False
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x_new = np.linspace(x.min(), x.max(), nx)
    y_new = np.linspace(y.min(), y.max(), ny)

    if nearest:
        x_idx = np.clip(np.searchsorted(x, x_new), 0, len(x) - 1)
        y_idx = np.clip(np.searchsorted(y, y_new), 0, len(y) - 1)
        z_new = z[np.ix_(y_idx, x_idx)]
        return x_new, y_new, z_new

    z_x = np.vstack([np.interp(x_new, x, row) for row in z])
    z_xy = np.vstack([np.interp(y_new, y, z_x[:, j]) for j in range(z_x.shape[1])]).T
    return x_new, y_new, z_xy


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
        "state": np.full(shape, np.nan, dtype=float),
    }

    limits = payload["limits"]
    for pt in filtered:
        row = int(np.where(flows == float(pt["totalFlow_Lpm"]))[0][0])
        col = int(np.where(voltages == float(pt["voltage_V"]))[0][0])
        grids["QtoBed_W"][row, col] = float(pt["QtoBed_W"])
        grids["totalCurrent_A"][row, col] = float(pt["totalCurrent_A"])
        grids["airOutlet_C"][row, col] = float(pt["airOutlet_C"])
        grids["wireMax_C"][row, col] = float(pt["wireMax_C"])
        grids["deltaP_Pa"][row, col] = float(pt["deltaP_Pa"])
        grids["state"][row, col] = feasibility_state_code(pt, limits)

    return grids, voltages, flows


def feasibility_state_code(point: dict, limits: dict) -> int:
    if bool(point["isFeasible"]):
        return 1
    if float(point["totalCurrent_A"]) > float(limits["maxTotalCurrent_A"]):
        return 2
    if float(point["airOutlet_C"]) > float(limits["maxAirOutletTemp_C"]):
        return 3
    if float(point["wireMax_C"]) > float(limits["maxWireTemp_C"]):
        return 4
    if float(point["deltaP_Pa"]) > float(limits["maxPressureDrop_Pa"]):
        return 5
    if float(point.get("waterLoss_kg_day", 0.0)) > float(limits.get("maxWaterLoss_kg_day", np.inf)):
        return 6
    if float(point.get("holeVelocity_m_s", 0.0)) > float(limits.get("maxHoleVelocity_m_s", np.inf)):
        return 7
    return 8


def plot_configuration_comparison(payload: dict, output_dir: Path) -> None:
    cfg = payload["configuration_comparisons"]
    labels = [item["label"] for item in cfg]
    colors = np.array([item["color"] for item in cfg], dtype=float)
    ua = np.array([float(item["UA_W_K"]) for item in cfg], dtype=float)
    heat = np.array([float(item["Qreq_W"]) for item in cfg], dtype=float)
    refs = []
    for item in cfg:
        rec = item.get("recommended", {})
        if isinstance(rec, dict) and rec.get("pitch_mm") is not None:
            refs.append(rec)
            continue
        refs.append(item.get("best_available", {}))
    power = np.array([float(item["totalPower_W"]) if item else np.nan for item in refs], dtype=float)
    flow = np.array([float(item["totalFlow_Lpm"]) if item else np.nan for item in refs], dtype=float)
    current = np.array([float(item["totalCurrent_A"]) if item else np.nan for item in refs], dtype=float)
    voltage = np.array([float(item["voltage_V"]) if item else np.nan for item in refs], dtype=float)
    pitch = np.array([float(item["pitch_mm"]) if item else np.nan for item in refs], dtype=float)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
    fig.suptitle(
        "Vessel Configuration Comparison\n"
        "Recommended = one comparison-selected operating point per configuration",
        fontsize=15,
        fontweight="bold",
    )

    x = np.arange(len(labels))
    panels = [
        (axes[0, 0], ua, "UA (W/K)", "Overall Ambient Conductance", "{:.2f}"),
        (axes[0, 1], heat, "Required heat (W)", "Heat Required at Setpoint", "{:.1f}"),
        (axes[1, 0], power, "Recommended power (W)", "Comparison-Based Recommended Power", "{:.1f}"),
    ]
    for ax, values, ylabel, title, fmt in panels:
        bars = ax.bar(x, values, color=colors)
        ax.set_xticks(x, labels, rotation=10)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        for bar, value in zip(bars, values):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height(),
                fmt.format(value),
                ha="center",
                va="bottom",
                fontsize=8,
            )

    ax = axes[1, 1]
    ax.bar(x, flow, color="#d4d4d8")
    ax.set_xticks(x, labels, rotation=10)
    ax.set_ylabel("Recommended flow (L/min)")
    ax.set_title("Flow, Current, Voltage, and Pitch")
    ax2 = ax.twinx()
    ax2.plot(x, current, "ko-", linewidth=1.5)
    ax2.set_ylabel("Recommended current (A)")
    for idx, (f_val, c_val, v_val, p_val) in enumerate(zip(flow, current, voltage, pitch)):
        ax.text(idx, f_val, f"{f_val:.0f} L/min\n{v_val:.0f} V, {p_val:.1f} mm", ha="center", va="bottom", fontsize=8)
        ax2.text(idx, c_val, f"{c_val:.2f} A", ha="center", va="bottom", fontsize=8)

    fig.text(
        0.02,
        0.02,
        payload["meta"]["recommended_definition"],
        fontsize=9,
    )
    fig.savefig(output_dir / "python_vessel_configuration_comparison.png", dpi=220)
    plt.close(fig)


def plot_constraint_maps(payload: dict, output_dir: Path, plot_points: int) -> None:
    grids, voltages, flows = build_active_grids(payload)
    rec = active_reference_point(payload)
    limits = payload["limits"]

    state_colors = np.array(
        [
            [0.10, 0.60, 0.20],
            [0.80, 0.20, 0.20],
            [0.93, 0.69, 0.13],
            [0.49, 0.18, 0.56],
            [0.00, 0.45, 0.74],
            [0.64, 0.08, 0.18],
            [0.85, 0.33, 0.10],
            [0.50, 0.50, 0.50],
        ]
    )
    state_cmap = matplotlib.colors.ListedColormap(state_colors)

    panels = [
        ("QtoBed_W", "Heat to Bed (W)", float(payload["requiredHeat_W"])),
        ("totalCurrent_A", "Total Current (A)", float(limits["maxTotalCurrent_A"])),
        ("airOutlet_C", "Outlet Air Temperature (C)", float(limits["maxAirOutletTemp_C"])),
        ("wireMax_C", "Max Wire Temperature (C)", float(limits["maxWireTemp_C"])),
        ("deltaP_Pa", "Total Pressure Drop (Pa)", float(limits["maxPressureDrop_Pa"])),
        ("state", "Constraint State Map", None),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
    fig.suptitle(
        f"Constraint and Performance Maps at {float(rec['pitch_mm']):.1f} mm Pitch\n"
        f"Star = comparison-based recommended point at {float(rec['voltage_V']):.1f} V and "
        f"{float(rec['totalFlow_Lpm']):.1f} L/min",
        fontsize=15,
        fontweight="bold",
    )

    for ax, (field, title, contour_value) in zip(axes.flat, panels):
        nearest = field == "state"
        x_new, y_new, z_new = resample_regular_grid(
            voltages,
            flows,
            grids[field],
            plot_points,
            plot_points,
            nearest=nearest,
        )
        extent = [x_new.min(), x_new.max(), y_new.min(), y_new.max()]
        if field == "state":
            im = ax.imshow(
                z_new,
                origin="lower",
                aspect="auto",
                extent=extent,
                cmap=state_cmap,
                vmin=1,
                vmax=8,
            )
            cb = fig.colorbar(im, ax=ax)
            cb.set_ticks([1, 2, 3, 4, 5, 6, 7, 8])
            cb.set_ticklabels(["feasible", "current", "air out", "wire temp", "pressure", "water loss", "hole vel", "heat shortfall"])
        else:
            im = ax.imshow(z_new, origin="lower", aspect="auto", extent=extent, cmap="viridis")
            fig.colorbar(im, ax=ax)
            if contour_value is not None:
                x_mesh, y_mesh = np.meshgrid(x_new, y_new)
                ax.contour(x_mesh, y_mesh, z_new, levels=[contour_value], colors="white", linewidths=1.0)

        ax.plot(float(rec["voltage_V"]), float(rec["totalFlow_Lpm"]), marker="*", markersize=14, color="white", markeredgecolor="black")
        ax.axvline(float(rec["voltage_V"]), color="white", linestyle=":", linewidth=1.0)
        ax.axhline(float(rec["totalFlow_Lpm"]), color="white", linestyle=":", linewidth=1.0)
        ax.set_xlabel("Voltage per branch (V)")
        ax.set_ylabel("Total flow (L/min)")
        ax.set_title(title)

    fig.text(
        0.02,
        0.02,
        payload["meta"]["star_definition"],
        fontsize=9,
    )
    fig.savefig(output_dir / "python_constraint_performance_maps_300.png", dpi=220)
    plt.close(fig)


def plot_wire_trade_study(payload: dict, output_dir: Path, plot_points: int) -> None:
    wire = payload["wire_diameter_sweep"]
    fig, axes = plt.subplots(3, 1, figsize=(14, 11), constrained_layout=True)
    fig.suptitle(
        "Wire Diameter Trade Study\n"
        "300-point interpolated display; targets come from vessel-comparison recommendations",
        fontsize=15,
        fontweight="bold",
    )

    series_defs = [
        ("optimalPower_W", "targetPower_W", "Matched total power (W)"),
        ("optimalFlow_Lpm", "targetFlow_Lpm", "Matched flow rate (L/min)"),
        ("optimalCurrent_A", "targetCurrent_A", "Matched current (A)"),
    ]

    for entry in wire:
        diam = as_float_array(entry["diameter_mm"])
        feasible = np.array(entry["isFeasible"], dtype=bool)
        valid_idx = np.where(feasible)[0]
        if valid_idx.size < 2:
            continue

        x_valid = diam[valid_idx]
        color = np.array(entry["color"], dtype=float)
        raw_best_index = entry["bestIndex"]
        if raw_best_index is None:
            best_index = None
        else:
            best_value = float(raw_best_index)
            best_index = int(best_value) - 1 if np.isfinite(best_value) else None
        x_dense = np.linspace(x_valid.min(), x_valid.max(), plot_points)

        for ax, (field, target_field, ylabel) in zip(axes, series_defs):
            y = as_float_array(entry[field])[valid_idx]
            interp = interp1d(x_valid, y, kind="linear", fill_value="extrapolate")
            y_dense = interp(x_dense)
            ax.plot(x_dense, y_dense, color=color, linewidth=1.8, label=entry["label"])
            ax.plot(x_valid, y, "o", color=color, markersize=4)
            target_value = entry[target_field]
            if target_value is not None:
                ax.axhline(float(target_value), color=color, linestyle="--", linewidth=1.0)
            if best_index is not None and best_index >= 0:
                ax.plot(
                    diam[best_index],
                    float(entry[field][best_index]),
                    marker="*",
                    markersize=12,
                    color="yellow",
                    markeredgecolor="black",
                )
            ax.set_ylabel(ylabel)
            ax.grid(True, alpha=0.3)

    axes[0].set_title("Power at Feasible Point Matched to Comparison Recommendation")
    axes[1].set_title("Flow at Feasible Point Matched to Comparison Recommendation")
    axes[2].set_title("Current at Feasible Point Matched to Comparison Recommendation")
    axes[2].set_xlabel("Wire diameter (mm)")
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        axes[0].legend(loc="best")
    fig.savefig(output_dir / "python_wire_diameter_trade_study_300.png", dpi=220)
    plt.close(fig)


def write_summary(payload: dict, output_dir: Path, plot_points: int, model_grid_points: int | None) -> None:
    rec = active_reference_point(payload)
    cfg = payload["active_configuration"]
    rec_kind = "Recommended point" if payload.get("recommended", {}).get("pitch_mm") is not None else "Best-available reference point"
    summary = textwrap.dedent(
        f"""
        Vermicomposter Heater Python Analysis
        ====================================

        Active configuration:
        - {cfg["label"]}

        {rec_kind} used throughout:
        - Pitch: {float(rec["pitch_mm"]):.1f} mm
        - Voltage: {float(rec["voltage_V"]):.1f} V
        - Total flow: {float(rec["totalFlow_Lpm"]):.1f} L/min
        - Total power: {float(rec["totalPower_W"]):.1f} W
        - Total current: {float(rec["totalCurrent_A"]):.2f} A

        What "recommended" means:
        - {payload["meta"]["recommended_definition"]}

        Why the same star appears on every constraint/performance plot:
        - {payload["meta"]["star_definition"]}
        - In this case the starred point is {float(rec["pitch_mm"]):.1f} mm, {float(rec["voltage_V"]):.1f} V, {float(rec["totalFlow_Lpm"]):.1f} L/min, so every panel shows the same x-y location while the z-value changes by metric.

        300-point display note:
        - Continuous Python plots were redrawn at {plot_points} display points per axis where that makes sense.
        - Vessel configuration bars remain discrete because there are only three actual vessel configurations in the model.
        - MATLAB model grid override: {model_grid_points if model_grid_points is not None else "not used"}.
        - {payload["meta"]["plot_note"]}
        """
    ).strip()

    (output_dir / "python_analysis_summary.txt").write_text(summary + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    output_dir = DEFAULT_OUTPUT_DIR.resolve()
    matlab_exe = find_matlab(args.matlab)
    study_config_path = args.study_config.resolve() if args.study_config is not None else None
    study_config = load_json(study_config_path)
    study_overrides = normalize_model_overrides_for_matlab(study_config.get("model_overrides", {}))
    matlab_parallel = matlab_parallel_config(study_config)
    export_json = export_results_json(
        matlab_exe, output_dir, args.model_grid_points, study_overrides, matlab_parallel
    )
    annotate_export_payload(export_json, study_overrides, study_config_path)
    synced_target = sync_export_to_study_data(export_json, study_config_path, study_config)
    payload = load_payload(export_json)

    plot_configuration_comparison(payload, output_dir)
    plot_constraint_maps(payload, output_dir, args.plot_points)
    plot_wire_trade_study(payload, output_dir, args.plot_points)
    write_summary(payload, output_dir, args.plot_points, args.model_grid_points)

    if synced_target is not None:
        print(f"Wrote analysis outputs to {output_dir} and synced export JSON to {synced_target}")
    else:
        print(f"Wrote analysis outputs to {output_dir}")


if __name__ == "__main__":
    main()

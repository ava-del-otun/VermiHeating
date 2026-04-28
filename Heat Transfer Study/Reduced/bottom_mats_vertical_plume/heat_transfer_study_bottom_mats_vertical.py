
#!/usr/bin/env python3
"""Bottom-mat heating + recreated cooling runner.

This script creates a new reduced model that:
1) keeps the recreated cooling workflow, and
2) replaces heating with four bottom mats plus a vertical porous-bed plume model.
3) inherits the recreated year-round Section 18 layer, including seeded
   spell-level heat-wave anomaly sampling.

Heating-side physics includes:
- buoyancy-driven superficial upward velocity closed by Ergun pressure-drop balance,
- passive stack-throughflow limit from bottom/top opening areas,
- Sherwood-based mass transfer for evaporation sink,
- 2D (length x height) transient-relaxation to steady temperature field.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


INCH_TO_M = 0.0254
R_V_WATER = 461.5  # J/(kg-K)
G_STD = 9.81  # m/s^2
MOLAR_MASS_O2_G_MOL = 32.0
SECONDS_PER_DAY = 86400.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Bottom-mat reduced model runner.")
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(__file__).with_name("study_config_bottom_mats_vertical.json"),
        help="Path to JSON config.",
    )
    return parser.parse_args()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def resolve_path(base_dir: Path, raw_path: str | Path) -> Path:
    p = Path(raw_path)
    return (base_dir / p).resolve() if not p.is_absolute() else p.resolve()


def json_default(obj: Any) -> Any:
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2, default=json_default), encoding="utf-8")


def load_base_module(base_script_path: Path):
    spec = importlib.util.spec_from_file_location("heat_transfer_study_recreated_base", str(base_script_path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load base script: {base_script_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def enabled_run_modes(config: dict) -> dict[str, bool]:
    run_modes = dict(config.get("run_modes", {}))
    return {
        "heating": bool(run_modes.get("heating", True)),
        "cooling": bool(run_modes.get("cooling", True)),
        "year_round": bool(run_modes.get("year_round", False)),
    }


def saturation_vapor_density_kg_m3(temp_c: np.ndarray | float) -> np.ndarray | float:
    """Tetens-based saturation vapor pressure converted to vapor density."""
    temp_c_arr = np.asarray(temp_c, dtype=float)
    p_sat_pa = 610.94 * np.exp((17.625 * temp_c_arr) / np.maximum(temp_c_arr + 243.04, 1e-9))
    rho_sat = p_sat_pa / (R_V_WATER * np.maximum(temp_c_arr + 273.15, 1e-9))
    if np.isscalar(temp_c):
        return float(rho_sat)
    return rho_sat


def bottom_mats_mode_config(config: dict) -> dict:
    raw = dict(config.get("heating_mode", {}))
    mats = dict(raw.get("mats", {}))
    bed = dict(raw.get("bed", {}))
    mass_transfer = dict(raw.get("mass_transfer", {}))
    boundary = dict(raw.get("boundary", {}))
    solver = dict(raw.get("solver", {}))
    sweep = dict(raw.get("sweep", {}))
    ventilation = dict(raw.get("ventilation", {}))
    bioheat = dict(raw.get("bioheat", {}))
    bio_temp = dict(bioheat.get("temperatureResponse", {}))
    enthalpy_map = dict(bioheat.get("enthalpy_kJ_per_molO2_byScenario", {}))

    mats.setdefault("count", 4)
    mats.setdefault("power_W_each", 17.5)
    mats.setdefault("supplyVoltage_V", 120.0)
    mats.setdefault("mat_length_in", 20.75)
    mats.setdefault("mat_width_in", 10.0)
    mats.setdefault("layout_rows", 2)
    mats.setdefault("layout_cols", 2)
    mats.setdefault("coupling_efficiency", 0.95)
    mats.setdefault("surface_deltaT_C", 2.0)

    bed.setdefault("length_in", 48.0)
    bed.setdefault("width_in", 24.0)
    bed.setdefault("height_in", 24.0)
    bed.setdefault("porosity", 0.80)
    bed.setdefault("particleDiameter_m", 0.006)
    bed.setdefault("kEff_W_mK", 0.45)
    bed.setdefault("bulkDensity_kg_m3", 650.0)
    bed.setdefault("cp_J_kgK", 3200.0)
    bed.setdefault("surfaceAreaDensity_m2_m3", 240.0)

    mass_transfer.setdefault("ambientRelativeHumidity", 0.80)
    mass_transfer.setdefault("vaporDiffusivity_m2_s", 2.4e-5)
    mass_transfer.setdefault("latentHeat_J_kg", 2.45e6)
    mass_transfer.setdefault("evaporationAvailabilityFactor", 0.85)
    mass_transfer.setdefault("dynamicBulkRHEnabled", True)
    mass_transfer.setdefault("closedLidBulkRHFloor", 0.97)
    mass_transfer.setdefault("maxEvaporation_kg_day", None)
    mass_transfer.setdefault("energyAdvectionFactor", 0.05)
    mass_transfer.setdefault("latentRecoveryFraction", 0.85)
    mass_transfer.setdefault("netMoistureLossFraction", 0.15)
    mass_transfer.setdefault("ventedLatentRecoveryFraction", 0.15)
    mass_transfer.setdefault("ventedNetMoistureLossFraction", 0.85)

    boundary.setdefault("sideU_W_m2K", 0.28)
    boundary.setdefault("lidU_W_m2K", 0.22)
    boundary.setdefault("bottomU_W_m2K", 0.20)

    solver.setdefault("nx", 140)
    solver.setdefault("nz", 90)
    solver.setdefault("iterations", 2600)
    solver.setdefault("dt_s", 0.25)
    solver.setdefault("convergence_tol_C", 1.0e-4)
    solver.setdefault("quickScale", 0.55)
    solver.setdefault("quickIterationsScale", 0.35)
    solver.setdefault("relaxation", 0.35)
    solver.setdefault("temperatureClipMax_C", 95.0)
    solver.setdefault("temperatureClipMarginBelowAmbient_C", 20.0)

    sweep.setdefault("duty_min", 0.05)
    sweep.setdefault("duty_max", 1.0)
    sweep.setdefault("duty_count", 24)

    ventilation.setdefault("enabled", True)
    ventilation.setdefault("bottomOpeningArea_m2", 0.0020)
    ventilation.setdefault("topOpeningArea_m2", 0.0020)
    ventilation.setdefault("bottomDischargeCoefficient", 0.62)
    ventilation.setdefault("topDischargeCoefficient", 0.62)
    ventilation.setdefault("mixingVelocityScale_m_s", 0.0015)

    bioheat.setdefault("enabled", True)
    bioheat.setdefault("scenario", "base")
    bioheat.setdefault("ourInitial_mgO2_gVS_day", 18.9)
    bioheat.setdefault("ourMature_mgO2_gVS_day", 2.5)
    bioheat.setdefault("decayK_day_inv", 0.0353)
    bioheat.setdefault("operatingDay", 30.0)
    bioheat.setdefault("vsFractionWet_kgVS_per_kgWet", 0.18)
    bioheat.setdefault("manualOUR_mgO2_gVS_day", None)
    bioheat.setdefault("manualQbio_W_m3", None)
    bioheat.setdefault("qbioClip_W_m3", [0.0, 500.0])

    enthalpy_map.setdefault("conservative", 313.6)  # Harper et al. 9.8 MJ/kg O2
    enthalpy_map.setdefault("base", 376.8)  # midpoint of literature range used here
    enthalpy_map.setdefault("aggressive", 440.0)  # de Guardia et al. fixed value
    bioheat["enthalpy_kJ_per_molO2_byScenario"] = enthalpy_map

    bio_temp.setdefault("enabled", True)
    bio_temp.setdefault("q10", 2.0)
    bio_temp.setdefault("reference_C", 22.0)
    bio_temp.setdefault("min_factor", 0.35)
    bio_temp.setdefault("max_factor", 1.65)
    bioheat["temperatureResponse"] = bio_temp

    cfg = {
        "method": "bottom_mats_vertical_plume",
        "ambientAir_C": float(raw.get("ambientAir_C", -15.0)),
        "designBed_C": float(raw.get("designBed_C", 22.0)),
        "mats": mats,
        "bed": bed,
        "mass_transfer": mass_transfer,
        "boundary": boundary,
        "solver": solver,
        "sweep": sweep,
        "ventilation": ventilation,
        "bioheat": bioheat,
    }
    return cfg


def heating_constraints_config(config: dict, heating_cfg: dict) -> dict:
    raw = dict(config.get("heating_constraints", {}))
    design_c = float(heating_cfg["designBed_C"])
    raw.setdefault("minBottomTemp_C", design_c)
    raw.setdefault("minTopTemp_C", design_c - 0.5)
    raw.setdefault("maxBottomTemp_C", design_c + 4.5)
    raw.setdefault("maxTopTemp_C", design_c + 4.5)
    raw.setdefault("maxSurfaceTemp_C", 90.0)
    raw.setdefault("maxWaterLoss_kg_day", 12.0)
    raw.setdefault("maxPressureDrop_Pa", 350.0)
    return raw


def duty_array(heating_cfg: dict) -> np.ndarray:
    sw = heating_cfg["sweep"]
    n = max(2, int(sw["duty_count"]))
    return np.linspace(float(sw["duty_min"]), float(sw["duty_max"]), n, dtype=float)


def bed_geometry_m(heating_cfg: dict) -> tuple[float, float, float]:
    bed = heating_cfg["bed"]
    return (
        float(bed["length_in"]) * INCH_TO_M,
        float(bed["width_in"]) * INCH_TO_M,
        float(bed["height_in"]) * INCH_TO_M,
    )


def mat_footprint_m(heating_cfg: dict) -> tuple[float, float]:
    mats = heating_cfg["mats"]
    length_per_mat = float(mats["mat_length_in"]) * INCH_TO_M
    width_per_mat = float(mats["mat_width_in"]) * INCH_TO_M
    total_length = max(1, int(mats["layout_cols"])) * length_per_mat
    total_width = max(1, int(mats["layout_rows"])) * width_per_mat
    return total_length, total_width


def _quick_solver_overrides(heating_cfg: dict) -> dict:
    cfg = json.loads(json.dumps(heating_cfg))
    solver = cfg["solver"]
    scale = max(0.2, float(solver.get("quickScale", 0.55)))
    iter_scale = max(0.1, float(solver.get("quickIterationsScale", 0.35)))
    solver["nx"] = max(48, int(round(scale * int(solver["nx"]))))
    solver["nz"] = max(32, int(round(scale * int(solver["nz"]))))
    solver["iterations"] = max(350, int(round(iter_scale * int(solver["iterations"]))))
    return cfg


def _bioheat_temperature_factor(temp_c: np.ndarray | float, bio_cfg: dict) -> np.ndarray | float:
    temp_arr = np.asarray(temp_c, dtype=float)
    temp_model = dict(bio_cfg.get("temperatureResponse", {}))
    if not bool(temp_model.get("enabled", True)):
        factor = np.ones_like(temp_arr)
    else:
        q10 = max(float(temp_model.get("q10", 2.0)), 1.0e-9)
        t_ref = float(temp_model.get("reference_C", 22.0))
        f_min = float(temp_model.get("min_factor", 0.35))
        f_max = float(temp_model.get("max_factor", 1.65))
        factor = np.power(q10, (temp_arr - t_ref) / 10.0)
        factor = np.clip(factor, min(f_min, f_max), max(f_min, f_max))
    if np.isscalar(temp_c):
        return float(factor)
    return factor


def _bioheat_base_terms(heating_cfg: dict) -> dict:
    bio = dict(heating_cfg.get("bioheat", {}))
    bed = dict(heating_cfg["bed"])
    enabled = bool(bio.get("enabled", True))
    scenario = str(bio.get("scenario", "base")).strip().lower() or "base"
    enthalpy_map = dict(bio.get("enthalpy_kJ_per_molO2_byScenario", {}))
    h_kj_mol = float(enthalpy_map.get(scenario, enthalpy_map.get("base", 376.8)))
    h_kj_mol = max(h_kj_mol, 0.0)

    if not enabled:
        return {
            "enabled": False,
            "scenario": scenario,
            "enthalpy_kJ_per_molO2": h_kj_mol,
            "our_mgO2_gVS_day": None,
            "qbio_W_per_kgVS": 0.0,
            "qbio_W_m3_base": 0.0,
            "vsFractionWet_kgVS_per_kgWet": float(np.clip(float(bio.get("vsFractionWet_kgVS_per_kgWet", 0.18)), 0.0, 1.0)),
            "operatingDay": float(max(float(bio.get("operatingDay", 30.0)), 0.0)),
            "decayK_day_inv": float(max(float(bio.get("decayK_day_inv", 0.0353)), 0.0)),
        }

    manual_qbio = bio.get("manualQbio_W_m3", None)
    qbio_w_m3: float
    our_day: float | None = None
    qbio_w_per_kg_vs = 0.0

    if manual_qbio is not None:
        qbio_w_m3 = float(max(float(manual_qbio), 0.0))
    else:
        manual_our = bio.get("manualOUR_mgO2_gVS_day", None)
        if manual_our is not None:
            our_day = float(max(float(manual_our), 0.0))
        else:
            our_initial = float(max(float(bio.get("ourInitial_mgO2_gVS_day", 18.9)), 0.0))
            our_mature = float(max(float(bio.get("ourMature_mgO2_gVS_day", 2.5)), 0.0))
            decay_k = float(max(float(bio.get("decayK_day_inv", 0.0353)), 0.0))
            operating_day = float(max(float(bio.get("operatingDay", 30.0)), 0.0))
            our_day = our_mature + (our_initial - our_mature) * float(np.exp(-decay_k * operating_day))

        qbio_w_per_kg_vs = our_day * (h_kj_mol / MOLAR_MASS_O2_G_MOL) * 1000.0 / SECONDS_PER_DAY
        vs_fraction = float(np.clip(float(bio.get("vsFractionWet_kgVS_per_kgWet", 0.18)), 0.0, 1.0))
        rho_bulk = float(max(float(bed.get("bulkDensity_kg_m3", 0.0)), 0.0))
        qbio_w_m3 = qbio_w_per_kg_vs * rho_bulk * vs_fraction

    clip_raw = bio.get("qbioClip_W_m3", [0.0, 500.0])
    if isinstance(clip_raw, (list, tuple)) and len(clip_raw) >= 2:
        qmin = float(clip_raw[0])
        qmax = float(clip_raw[1])
    else:
        qmin, qmax = 0.0, 500.0
    qbio_w_m3 = float(np.clip(qbio_w_m3, min(qmin, qmax), max(qmin, qmax)))

    return {
        "enabled": True,
        "scenario": scenario,
        "enthalpy_kJ_per_molO2": h_kj_mol,
        "our_mgO2_gVS_day": our_day,
        "qbio_W_per_kgVS": float(qbio_w_per_kg_vs),
        "qbio_W_m3_base": qbio_w_m3,
        "vsFractionWet_kgVS_per_kgWet": float(np.clip(float(bio.get("vsFractionWet_kgVS_per_kgWet", 0.18)), 0.0, 1.0)),
        "operatingDay": float(max(float(bio.get("operatingDay", 30.0)), 0.0)),
        "decayK_day_inv": float(max(float(bio.get("decayK_day_inv", 0.0353)), 0.0)),
    }


def _bioheat_reference_power_w(heating_cfg: dict, reference_temp_c: float) -> dict:
    base_terms = _bioheat_base_terms(heating_cfg)
    length_m, width_m, height_m = bed_geometry_m(heating_cfg)
    bed_volume_m3 = float(length_m * width_m * height_m)
    if not base_terms["enabled"]:
        return {
            "QbioReference_W": 0.0,
            "qbioReference_W_m3": 0.0,
            "temperatureFactor": 0.0,
            **base_terms,
        }
    temp_factor = float(_bioheat_temperature_factor(reference_temp_c, heating_cfg.get("bioheat", {})))
    qbio_ref_w_m3 = float(base_terms["qbio_W_m3_base"]) * temp_factor
    return {
        "QbioReference_W": qbio_ref_w_m3 * bed_volume_m3,
        "qbioReference_W_m3": qbio_ref_w_m3,
        "temperatureFactor": temp_factor,
        **base_terms,
    }


def _ergun_velocity_from_buoyancy(
    delta_t_k: np.ndarray,
    height_m: float,
    porosity: float,
    particle_d_m: float,
    rho_air: float,
    mu_air: float,
    ambient_c: float,
) -> np.ndarray:
    beta = 1.0 / max(ambient_c + 273.15, 1e-9)
    d_t = np.maximum(delta_t_k, 0.0)
    dp_buoy = rho_air * G_STD * beta * d_t * max(height_m, 1e-9)
    eps = max(min(porosity, 0.99), 0.05)
    dp = max(particle_d_m, 1e-9)
    a = 150.0 * mu_air * (1.0 - eps) ** 2 / max(eps**3 * dp**2, 1e-12)
    b = 1.75 * rho_air * (1.0 - eps) / max(eps**3 * dp, 1e-12)
    q = np.maximum(dp_buoy / max(height_m, 1e-9), 0.0)
    if b > 0.0:
        disc = np.maximum(a * a + 4.0 * b * q, 0.0)
        u = (-a + np.sqrt(disc)) / (2.0 * b)
    else:
        u = q / max(a, 1e-12)
    return np.maximum(np.nan_to_num(u, nan=0.0, posinf=0.0, neginf=0.0), 0.0)


def _stack_opening_superficial_velocity(
    delta_t_k: np.ndarray,
    height_m: float,
    ambient_c: float,
    footprint_area_m2: float,
    vent_cfg: dict,
) -> np.ndarray:
    if not bool(vent_cfg.get("enabled", False)):
        return np.full_like(np.asarray(delta_t_k, dtype=float), np.inf, dtype=float)
    a_bottom = max(float(vent_cfg.get("bottomOpeningArea_m2", 0.0)), 0.0)
    a_top = max(float(vent_cfg.get("topOpeningArea_m2", 0.0)), 0.0)
    cd_bottom = max(float(vent_cfg.get("bottomDischargeCoefficient", 0.62)), 1.0e-9)
    cd_top = max(float(vent_cfg.get("topDischargeCoefficient", 0.62)), 1.0e-9)
    if a_bottom <= 0.0 or a_top <= 0.0 or footprint_area_m2 <= 0.0:
        return np.zeros_like(np.asarray(delta_t_k, dtype=float), dtype=float)
    ae_bottom = cd_bottom * a_bottom
    ae_top = cd_top * a_top
    a_eff = 1.0 / np.sqrt(1.0 / max(ae_bottom**2, 1.0e-12) + 1.0 / max(ae_top**2, 1.0e-12))
    t_abs = max(float(ambient_c) + 273.15, 1.0)
    d_t = np.maximum(np.asarray(delta_t_k, dtype=float), 0.0)
    v_open = np.sqrt(np.maximum(2.0 * G_STD * max(height_m, 1.0e-9) * d_t / t_abs, 0.0))
    u_superficial = (a_eff / footprint_area_m2) * v_open
    return np.maximum(np.nan_to_num(u_superficial, nan=0.0, posinf=0.0, neginf=0.0), 0.0)

def simulate_bottom_mats_state(
    heating_cfg: dict,
    duty: float,
    ambient_c: float,
    quick: bool = False,
    include_field: bool = False,
) -> dict:
    cfg = _quick_solver_overrides(heating_cfg) if quick else heating_cfg
    mats = cfg["mats"]
    bed = cfg["bed"]
    mass = cfg["mass_transfer"]
    boundary = cfg["boundary"]
    solver = cfg["solver"]
    vent_cfg = cfg.get("ventilation", {})
    bio_cfg = cfg.get("bioheat", {})
    bio_terms = _bioheat_base_terms(cfg)
    qbio_base_w_m3 = float(bio_terms["qbio_W_m3_base"]) if bio_terms["enabled"] else 0.0

    length_m, width_m, height_m = bed_geometry_m(cfg)
    heat_len_m, heat_wid_m = mat_footprint_m(cfg)
    heat_len_m = min(heat_len_m, length_m)
    heat_wid_m = min(heat_wid_m, width_m)
    heat_area_m2 = max(heat_len_m * heat_wid_m, 1e-12)
    width_coverage = min(heat_wid_m / max(width_m, 1e-12), 1.0)
    bed_footprint_area_m2 = max(length_m * width_m, 1.0e-12)

    nx = max(24, int(solver["nx"]))
    nz = max(20, int(solver["nz"]))
    n_iter = max(120, int(solver["iterations"]))
    tol_c = max(float(solver["convergence_tol_C"]), 1.0e-8)
    relax = float(np.clip(float(solver.get("relaxation", 0.35)), 0.05, 1.0))
    temp_clip_max_c = float(solver.get("temperatureClipMax_C", 95.0))
    temp_clip_margin_below = max(float(solver.get("temperatureClipMarginBelowAmbient_C", 20.0)), 0.0)

    x = np.linspace(0.0, length_m, nx, dtype=float)
    z = np.linspace(0.0, height_m, nz, dtype=float)
    dx = length_m / max(nx - 1, 1)
    dz = height_m / max(nz - 1, 1)

    heat_start = 0.5 * (length_m - heat_len_m)
    heat_end = heat_start + heat_len_m
    heated_mask = (x >= heat_start) & (x <= heat_end)

    duty_clamped = float(np.clip(duty, 0.0, 1.0))
    p_total_w = (
        max(1, int(mats["count"]))
        * max(float(mats["power_W_each"]), 0.0)
        * duty_clamped
    )
    coupling = float(np.clip(float(mats["coupling_efficiency"]), 0.0, 1.0))
    q_to_bed_w = p_total_w * coupling
    q_flux_w_m2 = q_to_bed_w / heat_area_m2
    q_flux_x = np.where(heated_mask, q_flux_w_m2 * width_coverage, 0.0)

    k_eff = max(float(bed["kEff_W_mK"]), 1e-9)
    rho_bulk = max(float(bed["bulkDensity_kg_m3"]), 1e-9)
    cp_bulk = max(float(bed["cp_J_kgK"]), 1e-9)
    rho_cp = rho_bulk * cp_bulk
    alpha = k_eff / rho_cp
    d_particle = max(float(bed["particleDiameter_m"]), 1e-9)
    porosity = float(bed["porosity"])
    area_density = max(float(bed["surfaceAreaDensity_m2_m3"]), 0.0)

    # Air-side properties for section-12-style porous transport surrogates.
    t_ref_c = float(ambient_c)
    rho_air = 1.225 * 273.15 / max(t_ref_c + 273.15, 150.0)
    mu_air = 1.81e-5
    d_vapor = max(float(mass["vaporDiffusivity_m2_s"]), 1e-12)
    latent_j_kg = max(float(mass["latentHeat_J_kg"]), 1.0)
    rh_amb = float(np.clip(float(mass["ambientRelativeHumidity"]), 0.0, 1.0))
    dynamic_bulk_rh = bool(mass.get("dynamicBulkRHEnabled", True))
    rh_floor_closed = float(np.clip(float(mass.get("closedLidBulkRHFloor", 0.97)), 0.0, 1.0))
    max_evap_kg_day_raw = mass.get("maxEvaporation_kg_day", None)
    max_evap_kg_s = (
        max(float(max_evap_kg_day_raw), 0.0) / SECONDS_PER_DAY
        if max_evap_kg_day_raw not in (None, "")
        else None
    )
    energy_advection_factor = float(np.clip(float(mass.get("energyAdvectionFactor", 0.05)), 0.0, 1.0))
    vent_enabled = bool(vent_cfg.get("enabled", False))
    if vent_enabled:
        latent_recovery = float(np.clip(float(mass.get("ventedLatentRecoveryFraction", mass.get("latentRecoveryFraction", 0.15))), 0.0, 1.0))
        net_moisture_loss_fraction = float(
            np.clip(float(mass.get("ventedNetMoistureLossFraction", mass.get("netMoistureLossFraction", 0.85))), 0.0, 1.0)
        )
    else:
        latent_recovery = float(np.clip(float(mass.get("latentRecoveryFraction", 0.85)), 0.0, 1.0))
        net_moisture_loss_fraction = float(np.clip(float(mass.get("netMoistureLossFraction", 0.15)), 0.0, 1.0))
    vent_mix_velocity = max(float(vent_cfg.get("mixingVelocityScale_m_s", 0.0015)), 1.0e-9)
    rho_v_bulk_const = rh_amb * float(saturation_vapor_density_kg_m3(t_ref_c))
    rho_v_amb_col = np.full(nx, rho_v_bulk_const, dtype=float)
    sc = mu_air / max(rho_air * d_vapor, 1e-12)

    u_side = max(float(boundary["sideU_W_m2K"]), 0.0)
    u_lid = max(float(boundary["lidU_W_m2K"]), 0.0)
    u_bottom = max(float(boundary.get("bottomU_W_m2K", 0.0)), 0.0)
    u_front_back = max(float(boundary.get("frontBackU_W_m2K", u_side)), 0.0)
    bi_side = u_side * dx / max(k_eff, 1e-12)
    bi_lid = u_lid * dz / max(k_eff, 1e-12)
    bi_bottom = u_bottom * dz / max(k_eff, 1e-12)
    wall_sink_k_s = 2.0 * u_front_back / max(width_m * rho_cp, 1e-12)

    t = np.full((nz, nx), t_ref_c + 0.2, dtype=float)
    u_sup = np.zeros(nx, dtype=float)
    u_stack_cap = np.full(nx, np.inf, dtype=float)
    k_mass = np.zeros(nx, dtype=float)

    dx2 = max(dx * dx, 1e-12)
    dz2 = max(dz * dz, 1e-12)
    inv_laplacian_scale = (1.0 / dx2) + (1.0 / dz2)
    converged = False

    for it in range(n_iter):
        if it % 20 == 0:
            delta_col = np.maximum(t[0, :] - t[-1, :], 0.0)
            u_ergun = _ergun_velocity_from_buoyancy(
                delta_col,
                height_m,
                porosity,
                d_particle,
                rho_air,
                mu_air,
                t_ref_c,
            )
            u_stack_cap = _stack_opening_superficial_velocity(
                delta_col,
                height_m,
                t_ref_c,
                bed_footprint_area_m2,
                vent_cfg,
            )
            u_sup = np.minimum(u_ergun, u_stack_cap) if vent_enabled else u_ergun
            re = rho_air * np.maximum(u_sup, 1e-9) * d_particle / max(mu_air, 1e-12)
            sh = 2.0 + 1.1 * np.power(np.maximum(re, 1e-12), 0.6) * np.power(max(sc, 1e-12), 1.0 / 3.0)
            k_mass = sh * d_vapor / max(d_particle, 1e-12)

        rho_v_sat = saturation_vapor_density_kg_m3(t)
        if dynamic_bulk_rh:
            col_mean_c = np.asarray(np.mean(t, axis=0), dtype=float)
            rho_v_local_col = np.maximum(
                rho_v_amb_col,
                rh_floor_closed * np.asarray(saturation_vapor_density_kg_m3(col_mean_c), dtype=float),
            )
            if vent_enabled:
                vent_mix = np.clip(u_sup / np.maximum(u_sup + vent_mix_velocity, 1.0e-12), 0.0, 1.0)
                rho_v_bulk_col = (1.0 - vent_mix) * rho_v_local_col + vent_mix * rho_v_amb_col
            else:
                rho_v_bulk_col = rho_v_local_col
        else:
            rho_v_bulk_col = rho_v_amb_col
        evap_rate_vol_raw = k_mass[None, :] * area_density * np.maximum(rho_v_sat - rho_v_bulk_col[None, :], 0.0)
        lower_t = np.vstack([t[0:1, :], t[:-1, :]])
        grad_from_below = np.maximum((lower_t - t) / max(dz, 1e-12), 0.0)
        q_avail_w_m3 = rho_cp * np.maximum(u_sup[None, :], 0.0) * grad_from_below
        evap_availability = float(np.clip(float(mass.get("evaporationAvailabilityFactor", 0.85)), 0.0, 1.0))
        evap_w_m3 = np.minimum(latent_j_kg * evap_rate_vol_raw, evap_availability * q_avail_w_m3)
        if max_evap_kg_s is not None:
            evap_rate_total_kg_s_iter = float(np.sum(evap_w_m3 / max(latent_j_kg, 1e-12)) * dx * dz * width_m)
            if evap_rate_total_kg_s_iter > max_evap_kg_s > 0.0:
                evap_w_m3 *= max_evap_kg_s / evap_rate_total_kg_s_iter
        evap_w_m3_effective = (1.0 - latent_recovery) * evap_w_m3
        evap_rate_vol = evap_w_m3_effective / max(latent_j_kg, 1e-12)
        evap_sink_k_s = evap_w_m3_effective / max(rho_cp, 1e-12)
        t_new = np.array(t, copy=True)
        core = t[1:-1, 1:-1]
        u_up = np.maximum(energy_advection_factor * u_sup[None, 1:-1], 0.0)
        if qbio_base_w_m3 > 0.0:
            col_mean_c = np.asarray(np.mean(t, axis=0), dtype=float)
            bio_factor_col = np.asarray(_bioheat_temperature_factor(col_mean_c, bio_cfg), dtype=float)
            bio_source_k_s = (qbio_base_w_m3 * bio_factor_col[None, 1:-1]) / max(rho_cp, 1e-12)
        else:
            bio_source_k_s = 0.0
        lap_num = alpha * (
            (t[1:-1, 2:] + t[1:-1, :-2]) / dx2
            + (t[2:, 1:-1] + t[:-2, 1:-1]) / dz2
        )
        adv_num = u_up * t[:-2, 1:-1] / max(dz, 1.0e-12)
        denom = (
            2.0 * alpha * inv_laplacian_scale
            + u_up / max(dz, 1.0e-12)
            + wall_sink_k_s
        )
        candidate_core = (
            lap_num
            + adv_num
            + wall_sink_k_s * t_ref_c
            - evap_sink_k_s[1:-1, 1:-1]
            + bio_source_k_s
        ) / np.maximum(denom, 1.0e-12)
        t_new[1:-1, 1:-1] = core + relax * (candidate_core - core)

        # Bottom boundary with localized mat heat flux (Neumann).
        t_new[0, :] = (t_new[1, :] + q_flux_x * dz / max(k_eff, 1e-12) + bi_bottom * t_ref_c) / (1.0 + bi_bottom)

        # Top boundary conduction/convection to ambient; passive venting is handled in transport closure.
        t_new[-1, :] = (t_new[-2, :] + bi_lid * t_ref_c) / (1.0 + bi_lid)

        # End-wall losses (side boundaries in length-wise cut).
        t_new[:, 0] = (t_new[:, 1] + bi_side * t_ref_c) / (1.0 + bi_side)
        t_new[:, -1] = (t_new[:, -2] + bi_side * t_ref_c) / (1.0 + bi_side)

        t_new = np.clip(t_new, t_ref_c - temp_clip_margin_below, temp_clip_max_c)

        if it % 50 == 0:
            delta = float(np.max(np.abs(t_new - t)))
            if delta < tol_c:
                converged = True
                t = t_new
                break
        t = t_new

    # Final transport diagnostics
    delta_col = np.maximum(t[0, :] - t[-1, :], 0.0)
    u_ergun = _ergun_velocity_from_buoyancy(
        delta_col,
        height_m,
        porosity,
        d_particle,
        rho_air,
        mu_air,
        t_ref_c,
    )
    u_stack_cap = _stack_opening_superficial_velocity(
        delta_col,
        height_m,
        t_ref_c,
        bed_footprint_area_m2,
        vent_cfg,
    )
    u_sup = np.minimum(u_ergun, u_stack_cap) if vent_enabled else u_ergun
    re = rho_air * np.maximum(u_sup, 1e-9) * d_particle / max(mu_air, 1e-12)
    sh = 2.0 + 1.1 * np.power(np.maximum(re, 1e-12), 0.6) * np.power(max(sc, 1e-12), 1.0 / 3.0)
    k_mass = sh * d_vapor / max(d_particle, 1e-12)
    rho_v_sat = saturation_vapor_density_kg_m3(t)
    if dynamic_bulk_rh:
        col_mean_c = np.asarray(np.mean(t, axis=0), dtype=float)
        rho_v_local_col = np.maximum(
            rho_v_amb_col,
            rh_floor_closed * np.asarray(saturation_vapor_density_kg_m3(col_mean_c), dtype=float),
        )
        if vent_enabled:
            vent_mix = np.clip(u_sup / np.maximum(u_sup + vent_mix_velocity, 1.0e-12), 0.0, 1.0)
            rho_v_bulk_col = (1.0 - vent_mix) * rho_v_local_col + vent_mix * rho_v_amb_col
        else:
            rho_v_bulk_col = rho_v_local_col
    else:
        rho_v_bulk_col = rho_v_amb_col
    evap_rate_vol_raw = k_mass[None, :] * area_density * np.maximum(rho_v_sat - rho_v_bulk_col[None, :], 0.0)
    lower_t = np.vstack([t[0:1, :], t[:-1, :]])
    grad_from_below = np.maximum((lower_t - t) / max(dz, 1e-12), 0.0)
    q_avail_w_m3 = rho_cp * np.maximum(u_sup[None, :], 0.0) * grad_from_below
    evap_availability = float(np.clip(float(mass.get("evaporationAvailabilityFactor", 0.85)), 0.0, 1.0))
    evap_w_m3 = np.minimum(latent_j_kg * evap_rate_vol_raw, evap_availability * q_avail_w_m3)
    if max_evap_kg_s is not None:
        evap_rate_total_kg_s_iter = float(np.sum(evap_w_m3 / max(latent_j_kg, 1e-12)) * dx * dz * width_m)
        if evap_rate_total_kg_s_iter > max_evap_kg_s > 0.0:
            evap_w_m3 *= max_evap_kg_s / evap_rate_total_kg_s_iter
    evap_rate_vol_effective = ((1.0 - latent_recovery) * evap_w_m3) / max(latent_j_kg, 1e-12)
    evap_rate_total_kg_s = float(np.sum(evap_rate_vol_effective) * dx * dz * width_m)
    latent_w = evap_rate_total_kg_s * latent_j_kg
    evap_rate_vol_gross = evap_w_m3 / max(latent_j_kg, 1e-12)
    water_loss_total_kg_s = float(np.sum(evap_rate_vol_gross) * dx * dz * width_m * net_moisture_loss_fraction)

    if qbio_base_w_m3 > 0.0:
        col_mean_c = np.asarray(np.mean(t, axis=0), dtype=float)
        bio_temp_factor_col = np.asarray(_bioheat_temperature_factor(col_mean_c, bio_cfg), dtype=float)
        bio_temp_factor_field = np.repeat(bio_temp_factor_col[None, :], nz, axis=0)
        qbio_w_m3_field = qbio_base_w_m3 * bio_temp_factor_field
    else:
        bio_temp_factor_field = np.zeros_like(t)
        qbio_w_m3_field = np.zeros_like(t)
    qbio_total_w = float(np.sum(qbio_w_m3_field) * dx * dz * width_m)
    qbio_mean_w_m3 = float(np.mean(qbio_w_m3_field))
    qto_bed_total_w = float(q_to_bed_w + qbio_total_w)

    eps = max(min(porosity, 0.99), 0.05)
    a_erg = 150.0 * mu_air * (1.0 - eps) ** 2 / max(eps**3 * d_particle**2, 1e-12)
    b_erg = 1.75 * rho_air * (1.0 - eps) / max(eps**3 * d_particle, 1e-12)
    dp_col = height_m * (a_erg * u_sup + b_erg * u_sup**2)

    bottom_profile = t[0, :]
    top_profile = t[-1, :]
    t_bottom_mean_field = float(np.mean(bottom_profile))
    t_top_mean_field = float(np.mean(top_profile))
    t_lumped_mean = float(np.mean(t))
    t_bottom_mean = t_lumped_mean
    t_top_mean = t_lumped_mean
    t_bottom_max = float(np.max(bottom_profile))
    t_top_min = float(np.min(top_profile))
    vent_flow_m3_s = float(width_m * np.sum(np.maximum(u_sup, 0.0)) * dx)
    total_current_a = p_total_w / max(float(mats["supplyVoltage_V"]), 1e-9)
    our_day = bio_terms.get("our_mgO2_gVS_day", None)
    our_day_out = None if our_day is None else float(our_day)

    state = {
        "method": "bottom_mats_vertical_plume",
        "matDuty": duty_clamped,
        "totalPower_W": float(p_total_w),
        "totalCurrent_A": float(total_current_a),
        "QtoBed_W": float(qto_bed_total_w),
        "QmatToBed_W": float(q_to_bed_w),
        "Qbio_W": float(qbio_total_w),
        "qBioBase_W_m3": float(qbio_base_w_m3),
        "qBioMean_W_m3": float(qbio_mean_w_m3),
        "bioheatEnabled": bool(bio_terms["enabled"]),
        "bioheatScenario": str(bio_terms["scenario"]),
        "bioheatEnthalpy_kJ_per_molO2": float(bio_terms["enthalpy_kJ_per_molO2"]),
        "bioheatOUR_mgO2_gVS_day": our_day_out,
        "bioheat_W_per_kgVS": float(bio_terms["qbio_W_per_kgVS"]),
        "bioheatOperatingDay": float(bio_terms["operatingDay"]),
        "bioheatDecayK_day_inv": float(bio_terms["decayK_day_inv"]),
        "bioheatVsFractionWet_kgVS_per_kgWet": float(bio_terms["vsFractionWet_kgVS_per_kgWet"]),
        "QtoAmbient_W": float(max(p_total_w - q_to_bed_w, 0.0)),
        "airOutlet_C": float(t_top_mean),
        "wireMax_C": float(t_bottom_max + float(mats["surface_deltaT_C"])),
        "wireMean_C": float(t_bottom_mean + float(mats["surface_deltaT_C"])),
        "waterLoss_kg_day": float(max(water_loss_total_kg_s, 0.0) * 86400.0),
        "waterLossGross_kg_day": float(max(float(np.sum(evap_rate_vol_gross) * dx * dz * width_m), 0.0) * 86400.0),
        "netMoistureLossFraction": float(net_moisture_loss_fraction),
        "latentRecoveryFraction": float(latent_recovery),
        "latentEvap_W": float(max(latent_w, 0.0)),
        "deltaP_Pa": float(np.max(dp_col)),
        "ergunPressureDropMean_Pa": float(np.mean(dp_col)),
        "superficialVelocityMean_m_s": float(np.mean(u_sup)),
        "superficialVelocityPeak_m_s": float(np.max(u_sup)),
        "flowPerHole_Lpm": 0.0,
        "holeVelocity_m_s": float(np.max(u_sup)),
        "TbottomFullPower_C": t_bottom_mean,
        "TtopFullPower_C": t_top_mean,
        "TbedLumped_C": float(t_lumped_mean),
        "fieldBottomMean_C": float(t_bottom_mean_field),
        "fieldTopMean_C": float(t_top_mean_field),
        "bottomMin_C": float(np.min(bottom_profile)),
        "bottomMax_C": t_bottom_max,
        "topMin_C": t_top_min,
        "topMax_C": float(np.max(top_profile)),
        "heatingGradient_C": 0.0,
        "fieldGradient_C": float(t_bottom_mean_field - t_top_mean_field),
        "pressureSolveConverged": bool(converged),
        "convergedIterations": int(it + 1),
        "electricalTopologyLabel": "bottom_mats_parallel",
        "stringTubeCounts": [int(mats["count"])],
        "ambientAir_C": float(t_ref_c),
        "heatedLength_m": float(heat_len_m),
        "heatedWidth_m": float(heat_wid_m),
        "heatedZoneStart_m": float(heat_start),
        "heatedZoneEnd_m": float(heat_end),
        "passiveVentEnabled": bool(vent_enabled),
        "passiveVentFlow_Lpm": float(max(vent_flow_m3_s, 0.0) * 60000.0),
        "stackVelocityLimitMean_m_s": float(np.mean(u_stack_cap)) if np.all(np.isfinite(u_stack_cap)) else float("nan"),
    }
    if include_field:
        state["temperatureField_C"] = t
        state["x_m"] = x
        state["z_m"] = z
        state["qFluxBottom_W_m2"] = q_flux_x
        state["superficialVelocityProfile_m_s"] = u_sup
        state["ergunPressureDropProfile_Pa"] = dp_col
        state["bioheatField_W_m3"] = qbio_w_m3_field
        state["bioheatTempFactor"] = bio_temp_factor_field
    return state

def estimate_required_heating_w(heating_cfg: dict, ambient_c: float) -> dict:
    length_m, width_m, height_m = bed_geometry_m(heating_cfg)
    side_area = 2.0 * (length_m + width_m) * height_m
    top_area = length_m * width_m
    bottom_area = top_area
    boundary = heating_cfg["boundary"]
    ua = (
        max(float(boundary["sideU_W_m2K"]), 0.0) * side_area
        + max(float(boundary["lidU_W_m2K"]), 0.0) * top_area
        + max(float(boundary["bottomU_W_m2K"]), 0.0) * bottom_area
    )
    gross_w = max(ua * (float(heating_cfg["designBed_C"]) - float(ambient_c)), 0.0)
    bio_ref = _bioheat_reference_power_w(heating_cfg, float(heating_cfg["designBed_C"]))
    bio_credit_w = float(max(float(bio_ref["QbioReference_W"]), 0.0))
    net_external_w = max(gross_w - bio_credit_w, 0.0)
    return {
        "grossLoss_W": float(gross_w),
        "bioheatCredit_W": float(bio_credit_w),
        "requiredExternal_W": float(net_external_w),
        "bioheatReference": bio_ref,
    }


def constraints_for_state(state: dict, limits: dict, required_heat_w: float) -> tuple[bool, dict]:
    qmat_to_bed = float(state.get("QmatToBed_W", state.get("QtoBed_W", 0.0)))
    duty_needed = required_heat_w / max(qmat_to_bed, 1e-12) if qmat_to_bed > 0 else float("inf")
    residuals = {
        "bottom_shortfall_C": max(0.0, float(limits["minBottomTemp_C"]) - float(state["TbottomFullPower_C"])),
        "top_shortfall_C": max(0.0, float(limits["minTopTemp_C"]) - float(state["TtopFullPower_C"])),
        "bottom_over_C": max(0.0, float(state["TbottomFullPower_C"]) - float(limits["maxBottomTemp_C"])),
        "top_over_C": max(0.0, float(state["TtopFullPower_C"]) - float(limits["maxTopTemp_C"])),
        "surface_over_C": max(0.0, float(state["wireMax_C"]) - float(limits["maxSurfaceTemp_C"])),
        "water_over_kg_day": max(0.0, float(state["waterLoss_kg_day"]) - float(limits["maxWaterLoss_kg_day"])),
        "pressure_over_Pa": max(0.0, float(state["deltaP_Pa"]) - float(limits["maxPressureDrop_Pa"])),
        "duty_needed_over_1": max(0.0, duty_needed - 1.0),
    }
    feasible = all(v <= 1e-9 for v in residuals.values())
    residuals["dutyNeeded"] = float(duty_needed)
    return feasible, residuals


def build_bottom_mats_heating_payload(config: dict) -> tuple[dict, dict]:
    heating_cfg = bottom_mats_mode_config(config)
    limits = heating_constraints_config(config, heating_cfg)
    ambient_c = float(heating_cfg["ambientAir_C"])
    duties = duty_array(heating_cfg)
    required_heat = estimate_required_heating_w(heating_cfg, ambient_c)
    required_heat_w = float(required_heat["requiredExternal_W"])

    design_points: list[dict] = []
    feasible_points: list[dict] = []
    for duty in duties:
        pt = simulate_bottom_mats_state(heating_cfg, float(duty), ambient_c, quick=True, include_field=False)
        feasible, residuals = constraints_for_state(pt, limits, required_heat_w)
        pt["requiredHeat_W"] = float(required_heat_w)
        pt["requiredHeatGross_W"] = float(required_heat["grossLoss_W"])
        pt["bioheatCredit_W"] = float(required_heat["bioheatCredit_W"])
        pt["dutyNeeded"] = float(residuals["dutyNeeded"])
        pt["constraintResiduals"] = residuals
        pt["constraintState"] = "feasible" if feasible else "limited"
        design_points.append(pt)
        if feasible:
            feasible_points.append(pt)

    if feasible_points:
        recommended_seed = min(feasible_points, key=lambda item: (float(item["totalPower_W"]), float(item["matDuty"])))
    else:
        def infeasible_penalty(item: dict) -> tuple[float, float, float]:
            res = dict(item.get("constraintResiduals", {}))
            hard_penalty = float(
                max(0.0, float(res.get("bottom_shortfall_C", 0.0)))
                + max(0.0, float(res.get("top_shortfall_C", 0.0)))
                + max(0.0, float(res.get("bottom_over_C", 0.0)))
                + max(0.0, float(res.get("top_over_C", 0.0)))
                + max(0.0, float(res.get("surface_over_C", 0.0)))
                + max(0.0, float(res.get("water_over_kg_day", 0.0)))
                + max(0.0, float(res.get("pressure_over_Pa", 0.0)))
                + 50.0 * max(0.0, float(res.get("duty_needed_over_1", 0.0)))
            )
            return (hard_penalty, float(item.get("totalPower_W", 0.0)), float(item.get("matDuty", 0.0)))

        recommended_seed = min(design_points, key=infeasible_penalty)

    best_available = max(design_points, key=lambda item: float(item["QtoBed_W"]))

    detailed = simulate_bottom_mats_state(
        heating_cfg,
        float(recommended_seed["matDuty"]),
        ambient_c,
        quick=False,
        include_field=True,
    )
    feasible, residuals = constraints_for_state(detailed, limits, required_heat_w)
    detailed["requiredHeat_W"] = float(required_heat_w)
    detailed["requiredHeatGross_W"] = float(required_heat["grossLoss_W"])
    detailed["bioheatCredit_W"] = float(required_heat["bioheatCredit_W"])
    detailed["dutyNeeded"] = float(residuals["dutyNeeded"])
    detailed["constraintResiduals"] = residuals
    detailed["constraintState"] = "feasible" if feasible else "limited"

    recommended = {
        k: v
        for k, v in detailed.items()
        if k
        not in {
            "temperatureField_C",
            "x_m",
            "z_m",
            "qFluxBottom_W_m2",
            "superficialVelocityProfile_m_s",
            "ergunPressureDropProfile_Pa",
        }
    }

    payload = {
        "active_configuration": {"label": "Bottom mats + vertical porous plume (passive bottom-top venting)"},
        "selection_method": "minimum-power feasible duty using net external requirement (UA loss minus microbial source)",
        "design_points": design_points,
        "recommended": recommended,
        "best_available": best_available,
        "method": "bottom_mats_vertical_plume",
        "description": (
            "Heating uses four 20.75x10 in mats at the base (2x2 touching footprint), passive bottom/top openings, "
            "with buoyancy-driven upward porous transport closed by Ergun friction and stack-limited vent throughput, "
            "Sherwood-based "
            "mass-transfer latent sink, plus a literature-calibrated volumetric microbial heat source."
        ),
        "geometry": {
            "bedLength_in": 48.0,
            "bedWidth_in": 24.0,
            "bedHeight_in": 24.0,
            "heatedFootprint_in": [41.5, 20.0],
        },
        "limits": limits,
        "requiredHeat_W_designAmbient": float(required_heat_w),
        "requiredHeatGross_W_designAmbient": float(required_heat["grossLoss_W"]),
        "bioheatCredit_W_designReference": float(required_heat["bioheatCredit_W"]),
        "bioheat_reference": required_heat["bioheatReference"],
    }
    return payload, detailed


def select_heatmap_reference_duty(payload: dict, heating_cfg: dict) -> tuple[float, str]:
    points = list(payload.get("design_points", []))
    if points:
        target_c = float(heating_cfg.get("designBed_C", 22.0))

        def score(pt: dict) -> tuple[float, float, float]:
            tb = float(pt.get("TbottomFullPower_C", np.nan))
            tt = float(pt.get("TtopFullPower_C", np.nan))
            if not (np.isfinite(tb) and np.isfinite(tt)):
                return (float("inf"), float("inf"), float("inf"))
            mean_t = 0.5 * (tb + tt)
            spread_t = abs(tb - tt)
            return (abs(mean_t - target_c), spread_t, float(pt.get("matDuty", 1.0)))

        selected = min(points, key=score)
        return float(selected.get("matDuty", 1.0)), "closest-to-design-mean-temperature point"

    rec = payload.get("recommended", {}) if isinstance(payload.get("recommended", {}), dict) else {}
    if rec:
        return float(rec.get("matDuty", 1.0)), "recommended point fallback"
    return 1.0, "default duty fallback"


def plot_heating_sweep(payload: dict, output_dir: Path) -> None:
    points = list(payload.get("design_points", []))
    if not points:
        return
    points.sort(key=lambda p: float(p.get("matDuty", 0.0)))
    duty = np.array([float(p.get("matDuty", np.nan)) for p in points], dtype=float)
    qbed = np.array([float(p.get("QtoBed_W", np.nan)) for p in points], dtype=float)
    qmat = np.array([float(p.get("QmatToBed_W", np.nan)) for p in points], dtype=float)
    qbio = np.array([float(p.get("Qbio_W", 0.0)) for p in points], dtype=float)
    power = np.array([float(p.get("totalPower_W", np.nan)) for p in points], dtype=float)
    tbot = np.array([float(p.get("TbottomFullPower_C", np.nan)) for p in points], dtype=float)
    ttop = np.array([float(p.get("TtopFullPower_C", np.nan)) for p in points], dtype=float)
    water = np.array([float(p.get("waterLoss_kg_day", np.nan)) for p in points], dtype=float)
    dp = np.array([float(p.get("deltaP_Pa", np.nan)) for p in points], dtype=float)
    vel = np.array([float(p.get("superficialVelocityMean_m_s", np.nan)) for p in points], dtype=float)

    fig, ax = plt.subplots(2, 2, figsize=(11.0, 7.6), constrained_layout=True)
    ax[0, 0].plot(duty, qbed, label="Total bed heat (mat+bio) (W)", color="tab:red")
    ax[0, 0].plot(duty, qmat, label="Mat heat to bed (W)", color="tab:orange", linestyle="-.")
    ax[0, 0].plot(duty, qbio, label="Bioheat (W)", color="tab:green", linestyle=":")
    ax[0, 0].plot(duty, power, label="Electrical power (W)", color="tab:blue", linestyle="--")
    ax[0, 0].set_xlabel("Mat duty")
    ax[0, 0].set_ylabel("Power (W)")
    ax[0, 0].grid(alpha=0.25)
    ax[0, 0].legend()

    ax[0, 1].plot(duty, tbot, label="Bottom mean", color="tab:orange")
    ax[0, 1].plot(duty, ttop, label="Top mean", color="tab:green")
    ax[0, 1].set_xlabel("Mat duty")
    ax[0, 1].set_ylabel("Temperature (C)")
    ax[0, 1].grid(alpha=0.25)
    ax[0, 1].legend()

    ax[1, 0].plot(duty, water, color="tab:cyan", label="Water loss (kg/day)")
    ax[1, 0].set_xlabel("Mat duty")
    ax[1, 0].set_ylabel("Water loss (kg/day)")
    ax[1, 0].grid(alpha=0.25)
    ax[1, 0].legend()

    ax[1, 1].plot(duty, dp, color="tab:purple", label="Ergun drop (Pa)")
    ax[1, 1].set_xlabel("Mat duty")
    ax[1, 1].set_ylabel("Ergun pressure drop (Pa)", color="tab:purple")
    ax[1, 1].tick_params(axis="y", labelcolor="tab:purple")
    ax[1, 1].grid(alpha=0.25)
    ax_vel = ax[1, 1].twinx()
    ax_vel.plot(duty, vel, color="tab:brown", linestyle="--", label="Superficial velocity (m/s)")
    ax_vel.set_ylabel("Superficial velocity (m/s)", color="tab:brown")
    ax_vel.tick_params(axis="y", labelcolor="tab:brown")
    lines = ax[1, 1].get_lines() + ax_vel.get_lines()
    ax[1, 1].legend(lines, [ln.get_label() for ln in lines], loc="best")

    fig.suptitle("Bottom-Mat Heating Sweep (Vertical Porous Plume Model)")
    fig.savefig(output_dir / "heating_bottom_mats_sweep.png", dpi=220)
    plt.close(fig)

def plot_lengthwise_sideview_heatmap(detailed: dict, output_dir: Path) -> None:
    t = np.asarray(detailed["temperatureField_C"], dtype=float)
    x_m = np.asarray(detailed["x_m"], dtype=float)
    z_m = np.asarray(detailed["z_m"], dtype=float)
    heat_start_m = float(detailed["heatedZoneStart_m"])
    heat_end_m = float(detailed["heatedZoneEnd_m"])

    x_in = x_m / INCH_TO_M
    z_in = z_m / INCH_TO_M

    fig, ax = plt.subplots(figsize=(12.0, 4.6), constrained_layout=True)
    mesh = ax.imshow(
        t,
        origin="lower",
        aspect="auto",
        extent=[float(x_in[0]), float(x_in[-1]), float(z_in[0]), float(z_in[-1])],
        cmap="inferno",
        interpolation="bilinear",
    )
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label("Temperature (C)")

    ax.axvspan(
        heat_start_m / INCH_TO_M,
        heat_end_m / INCH_TO_M,
        ymin=0.0,
        ymax=0.03,
        color="#00c2ff",
        alpha=0.8,
        label="Heated base footprint",
    )
    ax.set_xlabel("Lengthwise position (in)")
    ax.set_ylabel("Height above base (in)")
    ax.set_title("Lengthwise Side-View Bed Temperature Map (Bottom Mats, Vertical Plume)")
    ax.legend(loc="upper right")

    fig.savefig(output_dir / "heating_sideview_lengthwise_heatmap.png", dpi=240)
    plt.close(fig)


def plot_lengthwise_sideview_heatmaps_with_without_bio(
    with_bio: dict,
    without_bio: dict,
    output_dir: Path,
) -> None:
    t_with = np.asarray(with_bio["temperatureField_C"], dtype=float)
    t_without = np.asarray(without_bio["temperatureField_C"], dtype=float)
    x_m = np.asarray(with_bio["x_m"], dtype=float)
    z_m = np.asarray(with_bio["z_m"], dtype=float)
    heat_start_m = float(with_bio["heatedZoneStart_m"])
    heat_end_m = float(with_bio["heatedZoneEnd_m"])

    x_in = x_m / INCH_TO_M
    z_in = z_m / INCH_TO_M
    extent = [float(x_in[0]), float(x_in[-1]), float(z_in[0]), float(z_in[-1])]

    common_min = float(min(np.min(t_with), np.min(t_without)))
    common_max = float(max(np.max(t_with), np.max(t_without)))
    delta = t_with - t_without
    dmax = float(np.max(np.abs(delta)))

    for title, arr, fname in (
        ("Lengthwise Side-View Temperature (With Microbial Heating)", t_with, "heating_sideview_lengthwise_heatmap_with_bio.png"),
        ("Lengthwise Side-View Temperature (Without Microbial Heating)", t_without, "heating_sideview_lengthwise_heatmap_without_bio.png"),
    ):
        fig, ax = plt.subplots(figsize=(12.0, 4.6), constrained_layout=True)
        mesh = ax.imshow(
            arr,
            origin="lower",
            aspect="auto",
            extent=extent,
            cmap="inferno",
            vmin=common_min,
            vmax=common_max,
            interpolation="bilinear",
        )
        cbar = fig.colorbar(mesh, ax=ax)
        cbar.set_label("Temperature (C)")
        ax.axvspan(
            heat_start_m / INCH_TO_M,
            heat_end_m / INCH_TO_M,
            ymin=0.0,
            ymax=0.03,
            color="#00c2ff",
            alpha=0.8,
            label="Heated base footprint",
        )
        ax.set_xlabel("Lengthwise position (in)")
        ax.set_ylabel("Height above base (in)")
        ax.set_title(title)
        ax.legend(loc="upper right")
        fig.savefig(output_dir / fname, dpi=240)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(12.0, 4.6), constrained_layout=True)
    mesh = ax.imshow(
        delta,
        origin="lower",
        aspect="auto",
        extent=extent,
        cmap="coolwarm",
        vmin=-dmax,
        vmax=dmax,
        interpolation="bilinear",
    )
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label("Delta T (with bio - without bio) (C)")
    ax.axvspan(
        heat_start_m / INCH_TO_M,
        heat_end_m / INCH_TO_M,
        ymin=0.0,
        ymax=0.03,
        color="#00c2ff",
        alpha=0.8,
        label="Heated base footprint",
    )
    ax.set_xlabel("Lengthwise position (in)")
    ax.set_ylabel("Height above base (in)")
    ax.set_title("Microbial-Heating Effect (Temperature Difference)")
    ax.legend(loc="upper right")
    fig.savefig(output_dir / "heating_sideview_lengthwise_heatmap_bio_delta.png", dpi=240)
    plt.close(fig)


def write_heating_summary(
    payload: dict,
    detailed: dict,
    output_dir: Path,
    heating_cfg: dict,
    heatmap_duty: float | None = None,
    heatmap_duty_note: str = "",
) -> None:
    rec = payload["recommended"]
    mats = heating_cfg["mats"]
    bio_ref = payload.get("bioheat_reference", {})
    bio_our = bio_ref.get("our_mgO2_gVS_day", None)
    bio_our_str = "n/a" if bio_our is None else f"{float(bio_our):.3f}"
    lines = [
        "============================================================",
        "BOTTOM-MAT HEATING MODEL SUMMARY",
        "============================================================",
        "Heating method                      : bottom_mats_vertical_plume",
        "Model basis                         : single-lumped bed control state + 2D diagnostic field + buoyancy/Ergun upward porous transport + stack-limited venting + Sherwood latent sink + microbial q_bio source",
        f"Bed geometry                         : 48.0 in (L) x 24.0 in (W) x 24.0 in (H)",
        f"Mats                                 : {int(mats['count'])} x 20.75 in x 10.00 in, touching 2x2 footprint",
        "Heated footprint                     : 41.50 in x 20.00 in (base plane)",
        f"Mat power rating                     : {float(mats['power_W_each']):.2f} W each at {float(mats['supplyVoltage_V']):.1f} VAC",
        f"Total electrical at 100% duty        : {float(mats['count']) * float(mats['power_W_each']):.2f} W",
        "",
        "Microbial q_bio model (literature-backed):",
        f"  Enabled                            : {bool(rec.get('bioheatEnabled', False))}",
        f"  Scenario                           : {str(rec.get('bioheatScenario', 'n/a'))}",
        f"  Enthalpy per mol O2                : {float(bio_ref.get('enthalpy_kJ_per_molO2', 0.0)):.2f} kJ/mol O2",
        f"  OUR at operating day               : {bio_our_str} mg O2/g VS/day",
        f"  q_bio base volumetric              : {float(bio_ref.get('qbio_W_m3_base', 0.0)):.3f} W/m3",
        f"  q_bio reference volumetric (Tset)  : {float(bio_ref.get('qbioReference_W_m3', 0.0)):.3f} W/m3",
        f"  q_bio reference total (Tset)       : {float(payload.get('bioheatCredit_W_designReference', 0.0)):.3f} W",
        "",
        "Design-ambient heating requirement decomposition:",
        f"  Gross envelope loss (UA*DeltaT)    : {float(payload.get('requiredHeatGross_W_designAmbient', 0.0)):.3f} W",
        f"  Minus microbial credit             : {float(payload.get('bioheatCredit_W_designReference', 0.0)):.3f} W",
        f"  External mat heat required         : {float(payload.get('requiredHeat_W_designAmbient', 0.0)):.3f} W",
        "",
        "Recommended heating operating point:",
        f"  Duty                               : {float(rec['matDuty']):.3f}",
        f"  Electrical power                   : {float(rec['totalPower_W']):.2f} W",
        f"  Mat heat to bed                    : {float(rec.get('QmatToBed_W', 0.0)):.2f} W",
        f"  Bioheat contribution               : {float(rec.get('Qbio_W', 0.0)):.2f} W",
        f"  Total heat to bed                  : {float(rec['QtoBed_W']):.2f} W",
        f"  Lumped bed temperature             : {float(rec.get('TbedLumped_C', rec['TbottomFullPower_C'])):.2f} C",
        f"  Legacy bottom/top aliases          : {float(rec['TbottomFullPower_C']):.2f} / {float(rec['TtopFullPower_C']):.2f} C",
        f"  Physical field bottom/top means    : {float(rec.get('fieldBottomMean_C', np.nan)):.2f} / {float(rec.get('fieldTopMean_C', np.nan)):.2f} C",
        f"  Physical field gradient            : {float(rec.get('fieldGradient_C', np.nan)):.2f} C",
        f"  Peak Ergun pressure drop           : {float(rec['deltaP_Pa']):.2f} Pa",
        f"  Mean superficial velocity          : {float(rec['superficialVelocityMean_m_s']):.5f} m/s",
        f"  Passive vent throughput            : {float(rec.get('passiveVentFlow_Lpm', np.nan)):.2f} L/min",
        f"  Water loss estimate (net)          : {float(rec['waterLoss_kg_day']):.4f} kg/day",
        f"  Water loss estimate (gross evap)   : {float(rec.get('waterLossGross_kg_day', np.nan)):.4f} kg/day",
        f"  Net moisture-loss fraction         : {float(rec.get('netMoistureLossFraction', np.nan)):.3f}",
        f"  Latent recovery fraction (effective): {float(rec.get('latentRecoveryFraction', np.nan)):.3f}",
        f"  Latent sink estimate               : {float(rec['latentEvap_W']):.4f} W",
        (
            f"  Heat-map reference duty            : {float(heatmap_duty):.3f} "
            f"({heatmap_duty_note})"
            if heatmap_duty is not None
            else "  Heat-map reference duty            : n/a"
        ),
        "",
        "Generated files:",
        "  - heating_sideview_lengthwise_heatmap.png",
        "  - heating_sideview_lengthwise_heatmap_with_bio.png",
        "  - heating_sideview_lengthwise_heatmap_without_bio.png",
        "  - heating_sideview_lengthwise_heatmap_bio_delta.png",
        "  - heating_bottom_mats_sweep.png",
        "  - heating_payload_bottom_mats.json",
        "  - heating_sideview_field_bottom_mats.npz",
        "============================================================",
    ]
    (output_dir / "study_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_cooling_with_base(base, config: dict, output_dir: Path) -> dict:
    payload = base.build_summer_cooling_payload(config)
    for stale in (
        "constraint_performance_maps.png",
        "total_current_curves.png",
        "moisture_perforation_maps.png",
        "wire_diameter_trade_study.png",
        "coil_diameter_trade_study.png",
        "wire_length_trade_study.png",
        "vessel_configuration_comparison.png",
        "summer_spot_cooler_curves.png",
        "summer_vertical_column_curves.png",
        "cooling_pareto_front.png",
        "summer_evaporative_cooling_curves.png",
        "radial_air_cooling_profile.png",
        "habitat_exclusion_cross_section.png",
    ):
        (output_dir / stale).unlink(missing_ok=True)
    base.plot_summer_cooling_curves(payload, config, output_dir)
    base.plot_cooling_pareto_front(payload, config, output_dir)
    base.plot_radial_air_cooling_profile(payload, config, output_dir)
    base.plot_habitat_exclusion_cross_section(payload, config, output_dir)
    try:
        base.write_cooling_summary(payload, config, output_dir)
    except Exception as exc:
        fallback_lines = [
            "============================================================",
            "COOLING SUMMARY (FALLBACK)",
            "============================================================",
            "Cooling branch completed, but base summary writer raised an exception.",
            f"Exception: {exc}",
            "Plots were still generated from the recreated cooling model.",
            "============================================================",
        ]
        (output_dir / "study_summary.txt").write_text("\n".join(fallback_lines) + "\n", encoding="utf-8")
    return payload



def run_year_round_with_base(
    base,
    config: dict,
    output_dir: Path,
    heating_cfg: dict,
    heating_payload: dict,
    cooling_payload: dict,
) -> dict:
    old_heating_sim = base.simulate_heating_reference_point

    def patched_heating_sim(config_local: dict, point: dict, ambient_air_c: float, bed_temp_c: float) -> dict:
        duty = float(point.get("matDuty", point.get("duty", 1.0)))
        state = simulate_bottom_mats_state(heating_cfg, duty, ambient_air_c, quick=True, include_field=False)
        state["wireMax_C"] = float(state.get("wireMax_C", ambient_air_c))
        state["airOutlet_C"] = float(state.get("airOutlet_C", ambient_air_c))
        return state

    try:
        base.simulate_heating_reference_point = patched_heating_sim
        yr_config = json.loads(json.dumps(config))
        yr_config.setdefault("parallel", {}).setdefault("python", {})["enabled"] = False
        yr_config.setdefault("parallel", {}).setdefault("python", {}).setdefault("loops", {})[
            "yearRoundDailySimulation"
        ] = False
        yr_annual = yr_config.setdefault("year_round_analysis", {})
        yr_climate = yr_annual.setdefault("climate_data", {})
        yr_heat = yr_climate.setdefault("heat_wave", {})
        yr_heat.setdefault("useDerivedSpellAnomalyBounds", True)
        yr_heat.setdefault(
            "derivedSpellAnomalyBoundsPath",
            "Heat Transfer Study/data/climate/noaa_heatwave_spell_anomaly_bounds_USW00014740_2016_2025.json",
        )
        yr_heat.setdefault("spellAnomalySamplingModel", "truncated_normal_between_monthly_min_max")
        yr_heat.setdefault("spellAnomalyRngSeed", 4143)
        annual_payload = base.build_year_round_payload(heating_payload, cooling_payload, yr_config)
        base.plot_year_round_energy_cost(annual_payload, output_dir)
        base.write_year_round_summary(annual_payload, output_dir)
        return annual_payload
    finally:
        base.simulate_heating_reference_point = old_heating_sim


def plot_combined_heating_cooling_dashboard(
    base,
    config: dict,
    heating_payload: dict,
    cooling_payload: dict,
    output_dir: Path,
) -> None:
    cool_points = list(cooling_payload.get("design_points", []))
    mat_points = list(heating_payload.get("design_points", []))
    if not cool_points or not mat_points:
        return

    cool_points.sort(key=lambda pt: float(pt.get("totalFlow_Lpm", np.nan)))
    mat_points.sort(key=lambda pt: float(pt.get("matDuty", np.nan)))

    cool_flow = np.array([float(pt.get("totalFlow_Lpm", np.nan)) for pt in cool_points], dtype=float)
    cool_tbed = np.array(
        [0.5 * (float(pt.get("TbottomFullPower_C", np.nan)) + float(pt.get("TtopFullPower_C", np.nan))) for pt in cool_points],
        dtype=float,
    )
    cool_cap = np.array([max(-float(pt.get("QtoBed_W", np.nan)), 0.0) for pt in cool_points], dtype=float)
    cool_tout = np.array([float(pt.get("airOutlet_C", np.nan)) for pt in cool_points], dtype=float)
    cool_water = np.array([float(pt.get("waterLoss_kg_day", np.nan)) for pt in cool_points], dtype=float)
    cool_power = np.array(
        [
            max(float(pt.get("assistBlowerPower_W", 0.0)), 0.0) + max(float(pt.get("spotCoolerPower_W", 0.0)), 0.0)
            for pt in cool_points
        ],
        dtype=float,
    )

    mat_tbed = np.array(
        [0.5 * (float(pt.get("TbottomFullPower_C", np.nan)) + float(pt.get("TtopFullPower_C", np.nan))) for pt in mat_points],
        dtype=float,
    )
    mat_power = np.array([float(pt.get("totalPower_W", np.nan)) for pt in mat_points], dtype=float)
    mat_flow = np.array([float(pt.get("totalFlow_Lpm", 0.0)) for pt in mat_points], dtype=float)
    mat_cap = np.array([max(float(pt.get("QtoBed_W", np.nan)), 0.0) for pt in mat_points], dtype=float)
    mat_tout = np.array([float(pt.get("airOutlet_C", np.nan)) for pt in mat_points], dtype=float)
    mat_water = np.array([float(pt.get("waterLoss_kg_day", np.nan)) for pt in mat_points], dtype=float)

    fig, ax = plt.subplots(2, 2, figsize=(13.5, 8.8), constrained_layout=True)
    cool_flow_order = np.argsort(np.nan_to_num(cool_flow, nan=np.inf))
    mat_flow_order = np.argsort(np.nan_to_num(mat_flow, nan=np.inf))

    ax[0, 0].plot(
        cool_flow[cool_flow_order],
        cool_tbed[cool_flow_order],
        color="tab:orange",
        linewidth=2.0,
        label="bed temperature (cooling)",
    )
    ax[0, 0].plot(
        mat_flow[mat_flow_order],
        mat_tbed[mat_flow_order],
        color="tab:green",
        linewidth=1.9,
        linestyle="--",
        label="bed temperature (heating)",
    )
    ax[0, 0].set_xlabel("Airflow (L/min)")
    ax[0, 0].set_ylabel("Bed temperature (C)")
    ax[0, 0].grid(alpha=0.25)
    ax[0, 0].legend()

    ax[0, 1].plot(
        cool_flow[cool_flow_order],
        cool_cap[cool_flow_order],
        color="tab:blue",
        linewidth=2.0,
        label="Cooling capacity",
    )
    ax[0, 1].plot(
        mat_flow[mat_flow_order],
        mat_cap[mat_flow_order],
        color="tab:orange",
        linewidth=1.8,
        linestyle="--",
        label="Heating capacity",
    )
    ax[0, 1].set_xlabel("Airflow (L/min)")
    ax[0, 1].set_ylabel("Capacity (W)")
    ax[0, 1].grid(alpha=0.25)
    ax_out = ax[0, 1].twinx()
    ax_out.plot(
        cool_flow[cool_flow_order],
        cool_tout[cool_flow_order],
        color="tab:red",
        linestyle="--",
        linewidth=1.8,
        label="outlet air temperature (cooling)",
    )
    ax_out.plot(
        mat_flow[mat_flow_order],
        mat_tout[mat_flow_order],
        color="tab:purple",
        linestyle="-.",
        linewidth=1.6,
        label="outlet air temperature (heating)",
    )
    ax_out.set_ylabel("Outlet air temperature (C)")
    ax_out.tick_params(axis="y", labelcolor="tab:red")
    lines = ax[0, 1].get_lines() + ax_out.get_lines()
    ax[0, 1].legend(lines, [ln.get_label() for ln in lines], loc="best")

    ax[1, 0].plot(mat_tbed, mat_power, color="tab:purple", linewidth=2.1, label="Heating power")
    cool_bed_order = np.argsort(np.nan_to_num(cool_tbed, nan=np.inf))
    ax[1, 0].plot(
        cool_tbed[cool_bed_order],
        cool_power[cool_bed_order],
        color="tab:blue",
        linewidth=1.9,
        linestyle="--",
        label="Cooling power",
    )
    ax[1, 0].set_xlabel("Bed temperature (C)")
    ax[1, 0].set_ylabel("Electrical power (W)")
    ax[1, 0].grid(alpha=0.25)
    ax[1, 0].legend()

    order = np.argsort(np.nan_to_num(cool_water, nan=np.inf))
    x_water_c = cool_water[order]
    y_power_c = cool_power[order]
    y_flow_c = cool_flow[order]
    ax[1, 1].plot(x_water_c, y_power_c, color="tab:blue", linewidth=2.0, label="Cooling electrical power")
    mat_order = np.argsort(np.nan_to_num(mat_water, nan=np.inf))
    x_water_h = mat_water[mat_order]
    y_power_h = mat_power[mat_order]
    y_flow_h = mat_flow[mat_order]
    if np.any(np.isfinite(x_water_h) & np.isfinite(y_power_h)):
        ax[1, 1].plot(
            x_water_h,
            y_power_h,
            color="tab:purple",
            linewidth=1.9,
            linestyle="--",
            label="Heating electrical power",
        )
    ax[1, 1].set_xlabel("Moisture loss (kg/day)")
    ax[1, 1].set_ylabel("Cooling electrical power (W)", color="tab:blue")
    ax[1, 1].tick_params(axis="y", labelcolor="tab:blue")
    ax[1, 1].grid(alpha=0.25)
    ax_flow = ax[1, 1].twinx()
    ax_flow.plot(
        x_water_c,
        y_flow_c,
        color="tab:orange",
        linestyle="--",
        linewidth=1.8,
        label="Fan flow rate (cooling)",
    )
    if np.any(np.isfinite(x_water_h) & np.isfinite(y_flow_h)):
        ax_flow.plot(
            x_water_h,
            y_flow_h,
            color="tab:green",
            linestyle="-.",
            linewidth=1.7,
            label="Fan flow rate (heating)",
        )
    ax_flow.set_ylabel("Fan flow rate (L/min)", color="tab:orange")
    ax_flow.tick_params(axis="y", labelcolor="tab:orange")
    lines = ax[1, 1].get_lines() + ax_flow.get_lines()
    ax[1, 1].legend(lines, [ln.get_label() for ln in lines], loc="best")

    fig.suptitle("Combined Heating/Cooling Dashboard (Bottom Mats + Cooling Column)")
    fig.savefig(output_dir / "combined_heating_cooling_dashboard.png", dpi=240)
    plt.close(fig)


def write_root_summary(root_output_dir: Path, mode_outputs: dict[str, Path]) -> None:
    lines = [
        "============================================================",
        "BOTTOM-MAT REDUCED MODEL RUN INDEX",
        "============================================================",
    ]
    for key in ("heating", "cooling", "year_round"):
        if key in mode_outputs:
            lines.append(f"{key} folder: {mode_outputs[key]}")
    lines.append(f"combined dashboard: {root_output_dir / 'combined_heating_cooling_dashboard.png'}")
    if not mode_outputs:
        lines.append("No run modes were enabled.")
    lines.append("============================================================")
    (root_output_dir / "study_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

def main() -> None:
    args = parse_args()
    config_path = args.config.resolve()
    config = load_json(config_path)

    root_output_dir = resolve_path(config_path.parent, config.get("output_dir", "outputs_bottom_mats_vertical_plume"))
    root_output_dir.mkdir(parents=True, exist_ok=True)
    run_modes = enabled_run_modes(config)
    heating_cfg = bottom_mats_mode_config(config)
    mode_outputs: dict[str, Path] = {}

    base_module = None
    if run_modes.get("cooling", False) or run_modes.get("year_round", False):
        base_script = resolve_path(config_path.parent, config.get("base_recreated_script", "..\\heat_transfer_study_recreated.py"))
        base_module = load_base_module(base_script)

    heating_payload = None
    cooling_payload = None

    if run_modes.get("heating", False):
        out_heating = root_output_dir / "heating"
        out_heating.mkdir(parents=True, exist_ok=True)
        heating_payload, detailed = build_bottom_mats_heating_payload(config)
        plot_heating_sweep(heating_payload, out_heating)
        duty_ref, duty_ref_note = select_heatmap_reference_duty(heating_payload, heating_cfg)
        heating_cfg_with_bio = json.loads(json.dumps(heating_cfg))
        heating_cfg_with_bio.setdefault("bioheat", {})["enabled"] = True
        heating_cfg_without_bio = json.loads(json.dumps(heating_cfg))
        heating_cfg_without_bio.setdefault("bioheat", {})["enabled"] = False
        detailed_with_bio = simulate_bottom_mats_state(
            heating_cfg_with_bio,
            duty_ref,
            float(heating_cfg_with_bio["ambientAir_C"]),
            quick=False,
            include_field=True,
        )
        detailed_without_bio = simulate_bottom_mats_state(
            heating_cfg_without_bio,
            duty_ref,
            float(heating_cfg_without_bio["ambientAir_C"]),
            quick=False,
            include_field=True,
        )
        plot_lengthwise_sideview_heatmap(detailed_with_bio, out_heating)
        plot_lengthwise_sideview_heatmaps_with_without_bio(detailed_with_bio, detailed_without_bio, out_heating)
        write_heating_summary(
            heating_payload,
            detailed,
            out_heating,
            heating_cfg,
            heatmap_duty=duty_ref,
            heatmap_duty_note=duty_ref_note,
        )
        np.savez_compressed(
            out_heating / "heating_sideview_field_bottom_mats.npz",
            temperatureField_C=np.asarray(detailed_with_bio["temperatureField_C"], dtype=float),
            x_m=np.asarray(detailed_with_bio["x_m"], dtype=float),
            z_m=np.asarray(detailed_with_bio["z_m"], dtype=float),
            qFluxBottom_W_m2=np.asarray(detailed_with_bio["qFluxBottom_W_m2"], dtype=float),
            superficialVelocityProfile_m_s=np.asarray(detailed_with_bio["superficialVelocityProfile_m_s"], dtype=float),
            ergunPressureDropProfile_Pa=np.asarray(detailed_with_bio["ergunPressureDropProfile_Pa"], dtype=float),
            bioheatField_W_m3=np.asarray(detailed_with_bio.get("bioheatField_W_m3", np.zeros_like(np.asarray(detailed_with_bio["temperatureField_C"], dtype=float))), dtype=float),
            bioheatTempFactor=np.asarray(detailed_with_bio.get("bioheatTempFactor", np.zeros_like(np.asarray(detailed_with_bio["temperatureField_C"], dtype=float))), dtype=float),
            temperatureField_with_bio_C=np.asarray(detailed_with_bio["temperatureField_C"], dtype=float),
            temperatureField_without_bio_C=np.asarray(detailed_without_bio["temperatureField_C"], dtype=float),
            bioheatField_with_bio_W_m3=np.asarray(detailed_with_bio.get("bioheatField_W_m3", np.zeros_like(np.asarray(detailed_with_bio["temperatureField_C"], dtype=float))), dtype=float),
            bioheatField_without_bio_W_m3=np.asarray(detailed_without_bio.get("bioheatField_W_m3", np.zeros_like(np.asarray(detailed_without_bio["temperatureField_C"], dtype=float))), dtype=float),
            temperatureField_recommendedDuty_C=np.asarray(detailed["temperatureField_C"], dtype=float),
            recommendedDuty=np.asarray([float(detailed.get("matDuty", np.nan))], dtype=float),
            duty_ref=np.asarray([duty_ref], dtype=float),
        )
        write_json(out_heating / "heating_payload_bottom_mats.json", heating_payload)
        mode_outputs["heating"] = out_heating

    if run_modes.get("cooling", False):
        if base_module is None:
            raise RuntimeError("Cooling requested but base recreated module could not be loaded.")
        out_cooling = root_output_dir / "cooling"
        out_cooling.mkdir(parents=True, exist_ok=True)
        cooling_payload = run_cooling_with_base(base_module, config, out_cooling)
        mode_outputs["cooling"] = out_cooling

    if base_module is not None and heating_payload is not None and cooling_payload is not None:
        plot_combined_heating_cooling_dashboard(base_module, config, heating_payload, cooling_payload, root_output_dir)

    if run_modes.get("year_round", False):
        if base_module is None:
            raise RuntimeError("Year-round requested but base recreated module could not be loaded.")
        if heating_payload is None or cooling_payload is None:
            raise RuntimeError("Year-round mode requires both heating and cooling payloads.")
        out_year = root_output_dir / "year_round"
        out_year.mkdir(parents=True, exist_ok=True)
        annual_payload = run_year_round_with_base(base_module, config, out_year, heating_cfg, heating_payload, cooling_payload)
        write_json(out_year / "year_round_payload_bottom_mats.json", annual_payload)
        mode_outputs["year_round"] = out_year

    write_root_summary(root_output_dir, mode_outputs)
    print(f"Wrote bottom-mat reduced-model outputs to {root_output_dir}")


if __name__ == "__main__":
    main()

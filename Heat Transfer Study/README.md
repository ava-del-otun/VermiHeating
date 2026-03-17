# Heat Transfer Study

This folder contains a standalone Python post-processing workflow for the vermicomposter heating study, the summer plenum-fed spot-cooler plus assist-blower study, and a year-round representative-climate energy/cost analysis layer built from the current heating and cooling reference points.

`heat_transfer_study.py` reads:

- `study_config.json` for the enabled run toggles under `run_modes`
- `study_config.json` for plot settings, wording, optimization thresholds, and MATLAB-export override inputs
- `data/model_export_grid100.json` for the exported model dataset

To Run:

- Dual-mode default:
  - Leave `run_modes.heating = true` and `run_modes.cooling = true`
  - Run: `python heat_transfer_study.py`
  - The script writes `outputs/heating/*`, `outputs/cooling/*`, `outputs/year_round/*`, and `outputs/study_summary.txt`
- Heating mode only:
  - Set `run_modes.heating = true` and `run_modes.cooling = false`
  - In the `Heating` directory run: `python ".\analyze_vermicomposter_heater.py" --model-grid-points 100 --plot-points 100 --output-dir ".\python_heater_analysis_moisture" --study-config ".\Heat Transfer Study\study_config.json"`
  - Next, in `Heat Transfer Study` run: `python heat_transfer_study.py`
- Summer cooling mode:
  - Set `run_modes.heating = false` and `run_modes.cooling = true`
  - Run only: `python heat_transfer_study.py`
  - The summer mode is built directly in Python and does not call the MATLAB heater sweep.
  - The current default architecture is: `Whynter ARC-1030WN -> insulated metal plenum -> Goorui GHBH 1D7 34 1R5 side channel blower -> splitter / aeration tubes`

Key config fields:

- `run_modes.heating`: enable/disable the MATLAB-backed heating post-processing branch
- `run_modes.cooling`: enable/disable the standalone summer plenum-fed spot-cooler branch
- `operating_mode`: legacy single-mode fallback retained for backward compatibility if `run_modes` is removed
- `plot_points`: interpolation grid density for contouring
- `constraint_curve_count`: number of contour curves per continuous constraint/performance panel
- `wire_plot_points`: interpolation density for the wire-diameter traces
- `auto_refresh_export`: automatic MATLAB-export refresh settings used by `heat_transfer_study.py`
- `matlab_parallel`: runtime-only MATLAB sweep acceleration settings; when enabled, the export uses `parfor` across independent pitch-flow rows while preserving the sequential voltage warm start within each row
- `python_parallel`: runtime-only Python process-pool settings; when enabled, the standalone script parallelizes the independent cooling flow-point sweep and the independent day solves in the year-round analysis
- `model_overrides`: shared study inputs passed into the MATLAB export, including sweep bounds, explicit operating limits, lexicographic thermal thresholds, wire-study flow resolution, perforation geometry, release fraction, moisture assumptions, the water-loss limit, and calibration multipliers
- `model_overrides.environment`: shared ambient and pressure inputs for the heating export (`greenhouseAir_C`, `pressure_Pa`, `referencePlateLength_m`)
- `model_overrides.bin`: active vessel geometry/top condition for the MATLAB export; the current default is the closed-top vented case and the widened-bin override now lives here as `width_m`
- `model_overrides.bin_wall`: editable vessel-wall material stack for the MATLAB export (`sheetThickness_mm`, `sheetK_W_mK`, `insulationThickness_mm`, `insulationK_W_mK`, `internalH_W_m2K`)
- `model_overrides.heater_tube`: editable base heater-tube geometry (`length_m`, `ID_mm`, `OD_mm`, `insulationThickness_mm`, `wallK_W_mK`, `insulationK_W_mK`, `externalH_W_m2K`, `segments`, `extraMinorK`)
- `model_overrides.heaterTube`: editable coil-geometry overrides such as `coilMeanD_m` and `coilSpan_m`
- `model_overrides.aeration`: editable aeration-network geometry and process inputs (`nParallelTubes`, `length_m`, `ID_mm`, `OD_mm`, optional `headerEnabled` / `headerID_mm` / `headerMinorK`, splitter inputs `splitterInletID_mm`, `branchConnectorID_mm`, `branchConnectorLength_m`, `splitterBodyK`, plus `wallK_W_mK`, `segments`, `bedH_W_m2K`, `releaseFraction`, `endMinorK`). `splitterOutletCount` is now derived automatically from `nParallelTubes`.
- `cooling_mode`: summer ambient temperature, Whynter supply target/capacity/COP, optional rated-spec derivation fields (`ratedTotalCapacity_BTU_h`, `ratedAirflow_CFM`, `ratedDehumidification_pints_day`, `deriveTargetFromRatedSpecs`), the `assist_blower` block, the assist-blower-limited airflow sweep, inlet RH, surface temperature, summer hard constraints, and cooling-mode lexicographic order
- `cooling_mode.assist_blower`: downstream pressure-source block for the plenum architecture (`model`, `ratedFlow_CFM`, `ratedPressure_inH2O`, `shutoffPressure_inH2O`, `motorPower_W`, `reference`); the current default basis is the 3-phase Goorui `GHBH 1D7 34 1R5`, and the code currently uses a linear fan-curve approximation between shutoff pressure and the implied free-air point when only the shutoff and rated-duty points are available
- `cost_analysis`: electricity-price inputs for operating-cost calculations; the current default uses the Connecticut residential average price from the U.S. EIA and can be overridden directly in `USD/kWh`
- `year_round_analysis`: representative Connecticut climate-year inputs for the annual temperature/cost layer, including NOAA 2016-2025 statewide monthly means and standard deviations plus NOAA daily-station freeze and heat-wave statistics
- `radial_profile`: standalone post-processing controls for the radial air-cooling / worm-safe-distance plot
- `tube_layout`: standalone post-processing controls for tube spacing, cross-section layout, and overlap-aware worm-habitat exclusion volume
- `optimization`: lexicographic selection criteria, including the spread limit, minimum bottom/top temperatures, top-candidate count, optimization-table count, and the ordered priority list
- `summary_inputs`: static spec/context values used to build the MATLAB-style text report

Editable geometry/material inputs:

- `model_overrides.bin_wall` is the editable source for the vessel steel/insulation stack.
- `model_overrides.environment` is the editable source for the shared heating ambient, pressure, and external-plate reference length.
- `model_overrides.heater_tube` is the editable source for heater-tube length, ID, OD, and heater insulation.
- `model_overrides.heaterTube` is the editable source for coil mean diameter and coil span.
- `model_overrides.aeration` is the editable source for aeration-tube count, dimensions, and splitter/header geometry.
- `model_overrides.aeration.headerEnabled = false` means the model assumes there is no separate large upstream header; pressure drop then includes only the splitter body, the branch-connector friction, and the connector-to-branch area change.
- `summary_inputs.bin_wall`, `summary_inputs.heater_tube`, and `summary_inputs.aeration` are now fallback spec text only.
- When the report is written, the corresponding `model_overrides.*` blocks override those fallback values before the spec block is printed.

Sweep bounds can be edited in `model_overrides.sweep` using either explicit arrays such as
`voltage_V` / `totalFlow_Lpm` / `wireDiameter_mm` / `coilMeanD_mm` / `wireLength_m` or the simpler
`*_min`, `*_max`, `*_count` fields. Explicit operating limits live in
`model_overrides.limits`, and the MATLAB-side lexicographic thermal thresholds live in
`model_overrides.optimization`. The active vessel mode now comes from
`model_overrides.bin.openTop`, `ventHoleCount`, and `ventHoleDiameter_m`.
The duty-based reporting criterion is also editable there as
`model_overrides.sweep.targetDuty`; the current default is `0.75` for the
`Minimum total flow at ... duty` line and the summary note about whether the
selected point can meet that duty target.

Additional helical-geometry studies:

- `model_overrides.sweep.coilMeanD_mm_*` drives the coil-mean-diameter trade study.
- `model_overrides.sweep.wireLength_m_*` drives the requested-wire-length trade study.
- The helical geometry is physically constrained so the coil outer diameter cannot exceed the heater-tube inner diameter, i.e. `D_coil,mean + d_wire <= ID_heater`.
- The wire-length study is implemented as a helical sensitivity study: for each requested wire length, the code solves the implied helix pitch for the configured coil span and mean diameter rather than treating wire length as a purely arbitrary electrical override.

Distribution-manifold pressure model:

- The heating MATLAB model and the standalone cooling model now treat the common distribution hardware explicitly instead of using only a single lumped header `K`.
- The implemented path is:
  - optional upstream header loss,
  - optional header-to-splitter contraction,
  - splitter-body minor loss,
  - branch-connector friction,
  - connector-to-downstream contraction or sudden expansion.
- Sudden expansion uses Idel'chik Section IV sudden-expansion loss, written on the upstream-velocity basis as `zeta = (1 - A1/A2)^2`; this is algebraically the same Borda-Carnot coefficient used in White and Munson.
- Contraction uses Idel'chik Section III, Diagram 3-5, ``Conical converging bellmouth without end wall,'' evaluated with the editable cone angle and the relation `zeta = C (1/epsilon - 1)` on the smaller-pipe velocity basis.
- The current JSON default matches the stated hardware assumption:
  - no separate header,
  - `4` splitter outlets,
  - `19 mm` splitter inlet and branch connector,
  - `60 deg` conical contraction angle,
  - configurable connector length and splitter-body `K`.
- The GroundWork 4-way brass manifold basis is the Tractor Supply `4510201` product, which is a `3/4 in` 4-way shutoff manifold. The model converts that nominal `3/4 in` size to an editable `19 mm` hydraulic placeholder because the product page does not publish an internal hydraulic diameter or internal loss coefficient. The unresolved internal manifold loss is represented by `splitterBodyK`.

Automatic export refresh:

- If `auto_refresh_export.enabled` is true, `heat_transfer_study.py` compares the current
  `model_overrides` block against the signature stored in the export JSON.
- If the signatures differ, it automatically calls `analyze_vermicomposter_heater.py`
  before generating plots and the `.txt` summary.
- If `matlab_parallel.enabled` is true and MATLAB Parallel Computing Toolbox is available, the MATLAB export uses `parfor` with the configured worker count (`0` = local-profile default). This changes runtime behavior only; it does not change the physics or the export signature.
- If `python_parallel.enabled` is true, `heat_transfer_study.py` uses a Windows-safe process pool for the coarse standalone tasks that are actually independent: the summer cooling flow sweep and the representative-year daily solves. This also changes runtime behavior only; it does not change the physics or the export signature.
- The one exception is a pure vessel-top switch handled through the already-exported
  `configuration_comparisons` data: with
  `auto_refresh_export.allow_bin_only_switch_from_export = true`, the standalone study can
  switch between `open top`, `closed top`, and `closed top with vents` without rerunning MATLAB
  if those comparison datasets are already present in the export JSON and no other physics-side
  `model_overrides` entries changed.

Wetted-area model:

- By default the evaporation model derives wetted area from the particle-size distribution
  reported in Frederickson et al. (2007), Table 3.
- The default basis is `model_overrides.evaporation.wettedAreaBasis = "table3_particle_surface"`.
- In that mode the code blends the windrow and vermicomposted size fractions, treats each bin
  as an equivalent-sphere class with representative diameters from
  `model_overrides.evaporation.table3ParticleSurface.representativeDiameters_mm`, computes the
  Sauter mean diameter
  `d32 = 1 / sum(w_i / d_i)`, and then uses `A_wet = f_access * 6 * V_fill / d32`.
- `accessibleSurfaceFraction` is the exposed-wetted-surface fraction applied to that particle-surface estimate.
- Use `model_overrides.evaporation.wettedAreaBasis = "top_plan_area"` with
  `wettedAreaFraction` only if you want the older plan-area approximation.
- Use `wettedAreaBasis = "manual"` and provide `wettedArea_m2` only if you want a direct override.

Evaporation coupling:

- Evaporation is not just a post-processing plot. It changes the base aeration heat-transfer result through `Q_to_bed = sum(q_wall + q_jet - q_evap)`.
- The configured `aeration.bedH_W_m2K` remains an effective sensible wall-transfer coefficient.
- Latent evaporation is modeled separately through the moisture calculation and should not also be added into `bedH_W_m2K`.

Editing guidance:

- If `auto_refresh_export.enabled` is true, `heat_transfer_study.py` will rerun `analyze_vermicomposter_heater.py` automatically whenever `model_overrides.*` changes, except for the pure top-mode switch described above.
- The current default auto-refresh profile uses a reduced runtime setting: `60 x 60` for the main MATLAB export grid, reduced geometry-study sampling, and reduced heater/aeration segmentation.
- If `auto_refresh_export.enabled` is false, you must rerun `analyze_vermicomposter_heater.py` manually every time `model_overrides.*` changes.
- The one exception is the summer cooling branch enabled through `run_modes.cooling`: that mode is generated directly in Python, so no MATLAB refresh is needed.
- If you change only top-level `optimization`, `radial_profile`, wording fields, or plot-density fields, rerunning `heat_transfer_study.py` is enough.
- The top-level `optimization.priority_order` controls the standalone Python lexicographic ordering and report tables.
- The default clean ordering is: hard constraints first, then maximize `min(T_bottom, T_top)`, then minimize total power, then minimize spread, with mean bed temperature and heat shortfall used only as late tie-breakers.
- The MATLAB-side thermal thresholds used in the exported criteria logic live in `model_overrides.optimization`.

Radial air-cooling plot:

- `radial_profile.temperature_threshold_C` is the worm-safety threshold used to mark the safe distance.
- `radial_profile.initial_plume_radius_mode`, `initial_plume_radius_m`, `spread_rate_m_per_m`, and `max_plume_radius_m`
  control the finite plume-spreading model.
- The plot starts at the aeration-tube outer surface and applies a finite-radius top-hat plume model
  to the active recommended or best-available point.
- The plume energy balance is re-solved locally as
  `m'_rel c_p dT/dx = -(q'_sens + q'_lat)` with
  `q'_sens = h_sens pi a_wet b(x)^2 (T - T_sink)` and
  `q'_lat = m'_evap h_fg`.
- The local evaporation term is recomputed along the radial profile using the same IAPWS saturation properties and Lewis-relation mass-transfer framework as the main moisture model, so the latent sink falls as the plume cools and humidifies.

Tube layout and habitat exclusion:

- The default `tube_layout.layoutMode = "rule_based_1to12"` applies the requested orientation rules for 1 through 12 tubes inside the filled-bed rectangle, using 10% of the filled-bed width/height as the wall offset.
- For `8` through `12` tubes, the implemented rule is:
  - `8 = 3` on each side + `2` on the bottom
  - `9 = 3` on each side + `3` on the bottom
  - `10 = 3` on each side + `4` on the bottom
  - `11 = 4` on each side + `3` on the bottom
  - `12 = 4` on each side + `4` on the bottom
- `tube_layout` precedence is now explicit and consistent:
  - `tubeCenters_m` overrides everything
  - otherwise `tubeCenterX_m` uses a single row at `tubeCenterHeight_m`
  - otherwise `layoutMode` is used
- If `tubeCenterHeight_m` is left `null`, the code defaults it to the same `wallOffsetHeightFraction * fill_height` used by the rule-based bottom row.
- `model_overrides.aeration.nParallelTubes` controls how many tubes are used in both the MATLAB model and the rule-based layout.
- If `tube_layout.tubeCenters_m` is provided, those explicit `[x_m, y_m]` coordinates override the rule-based pattern.
- The habitat exclusion envelope uses the tube radius plus the plume-model safe-distance cutoff.
- The reported unsafe volume is the overlap-aware union of those exclusion cylinders plus the optional 1D wall-return bands from the left, right, and bottom vessel sheets, clipped to the filled-bed rectangle in cross section and then extruded over tube length.
- `tube_layout.cross_section_grid_points` controls the midpoint-grid resolution used for the overlap calculation.

Year-round climate inputs:

- `year_round_analysis.climate_data.monthlyMean_C` and `monthlyStd_C` are the Connecticut statewide 2016-2025 monthly mean-temperature and monthly standard-deviation arrays derived from NOAA NCEI Climate at a Glance.
- `year_round_analysis.climate_data.freeze.*` and `heat_wave.*` are derived from NOAA daily-summaries station `USW00014740` over the same 2016-2025 period.
- The current annual driver builds a representative 365-day year by:
  - reconstructing a smooth daily ambient baseline from the 12 monthly Connecticut mean-temperature anchors using periodic Gaussian kernel regression on day of year,
  - inserting freeze days with Gaussian timing weights centered on the observed mean freeze day and a repulsion penalty so the selected days are dispersed rather than packed into one block,
  - inserting heat-wave spells with Gaussian timing weights on spell starts plus repulsion between starts, while still keeping each individual heat-wave spell contiguous by definition,
  - then re-solving the selected heating and cooling operating points at each day's ambient temperature before computing daily duty, electricity, water use, and cost.
- The annual plot in `outputs\year_round\year_round_energy_cost_analysis.png` is now a time-series figure:
  - top panel: ambient, bottom-node, and top-node temperatures vs day of year
  - bottom panel: daily cost and cumulative cost vs day of year

Run it from `Heating` with:

```powershell
python ".\Heat Transfer Study\heat_transfer_study.py"
```

Outputs are written to `Heat Transfer Study\outputs`.

Root output set:

- `study_summary.txt` as an index pointing to the enabled mode-specific folders
- `heating\...` if `run_modes.heating = true`
- `cooling\...` if `run_modes.cooling = true`
- `year_round\...` if both heating and cooling are enabled and `year_round_analysis.enabled = true`

Heating-mode output set (`outputs\heating`):

- `constraint_performance_maps.png`
- `total_current_curves.png`
- `wire_diameter_trade_study.png`
- `coil_diameter_trade_study.png`
- `wire_length_trade_study.png`
- `moisture_perforation_maps.png`
- `radial_air_cooling_profile.png`
- `habitat_exclusion_cross_section.png`
- `study_summary.txt`

Summer-mode output set (`outputs\cooling`):

- `summer_spot_cooler_curves.png`
- `radial_air_cooling_profile.png`
- `habitat_exclusion_cross_section.png`
- `study_summary.txt`

Year-round output set (`outputs\year_round`):

- `year_round_energy_cost_analysis.png`
  - representative-year temperature and cost vs time
- `study_summary.txt`
  - source URLs, annual totals, monthly climate/event table, and monthly operating table

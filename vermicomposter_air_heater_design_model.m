function results = vermicomposter_air_heater_design_model(varargin)
% VERMICOMPOSTER AIR HEATER DESIGN MODEL
%
% Correlation/source map used in comments below
%   [R1] Churchill, S.W.; Bernstein, M. "A Correlating Equation for Forced
%        Convection from Gases and Liquids to a Circular Cylinder in
%        Crossflow." Journal of Heat Transfer 99(2), 1977, 300-306.
%   [R2] Gnielinski, V. "New Equations for Heat and Mass Transfer in
%        Turbulent Pipe and Channel Flow." Int. Chem. Eng. 16(2), 1976.
%   [R3] Churchill, S.W. "Friction-Factor Equation Spans All Fluid-Flow
%        Regimes." Chemical Engineering 84(24), 1977.
%   [R4] Churchill, S.W.; Chu, H.H.S. "Correlating Equations for Laminar
%        and Turbulent Free Convection from a Vertical Plate."
%        Int. J. Heat Mass Transfer 18(11), 1975.
%   [R5] Bergman, Lavine, Incropera, DeWitt. Fundamentals of Heat and Mass
%        Transfer. Thermal resistance and lumped-capacitance balances.
%   [R6] Idel'chik, I.E. Handbook of Hydraulic Resistance, Section IV
%        sudden-expansion loss. This is the same Borda-Carnot coefficient
%        used in White/Munson when referenced to the upstream velocity.
%   [R7] Ideal-gas density and Sutherland-law air viscosity constants.
%   [R8] Hartenstein, R. "Metabolic parameters of the earthworm Eisenia
%        foetida in relation to temperature." Biotechnol. Bioeng. 24(8),
%        1982.
%   [R9] IAPWS SR1-86(1992). Revised Supplementary Release on Saturation
%        Properties of Ordinary Water Substance.
%   [R10] Hasan, A.; Siren, K. "The Lewis factor and its influence on the
%         performance prediction of wet-cooling towers." International
%         Journal of Thermal Sciences 44(9), 2005.
%   [R11] Frederickson, J.; Howell, G.; Hobson, A.M. "Effect of
%         pre-composting and vermicomposting on compost characteristics."
%         European Journal of Soil Biology 43 (2007) S320-S326. Table 3 is
%         used here for the particle-size distribution of the matured
%         material.
%   [R12] Idel'chik, I.E. Handbook of Hydraulic Resistance. Sudden
%         contraction/expansion and junction-loss methods used here for the
%         splitter/connector distribution network.
%
% Notes
%   - The default maximum heater outlet temperature is reduced to 30 C
%     because local warm jets from the perforations can expose worms to
%     temperatures above the bulk-bed average. [R8]
%   - Adding parallel tubes improves distribution, but it does not change
%     the thermodynamic minimum flow implied by Q = m_dot * cp * DeltaT.

params = defaultParams();

%% USER INPUTS - edit these first
params.environment.greenhouseAir_C = -10;
params.bin.setpoint_C = 20;
params.bin.minSafe_C = 15;
params.bin.maxSafe_C = 25;

params.bin.length_m = 1.22;
params.bin.width_m = 0.61;
params.bin.height_m = 0.61;
params.bin.openTop = true;
params.bin.ventHoleCount = 1;
params.bin.ventHoleDiameter_m = 0.0;
params.bin.comparisonVentHoleDiameter_m = 0.025;
% For a covered top with a small vent, set:
%   params.bin.openTop = false;
%   params.bin.ventHoleDiameter_m = desired hole diameter (m);

params.heaterTube.length_m = 0.30;
params.heaterTube.ID_m = 0.015;
params.heaterTube.OD_m = 0.019;
params.heaterTube.coilMeanD_m = 0.0115;
params.heaterTube.pitch_m = 0.008;
params.heaterTube.coilSpan_m = params.heaterTube.length_m;
params.heaterTube.manualWireLength_m = NaN;

params.aeration.ID_m = 0.040;
params.aeration.OD_m = 0.044;
params.aeration.length_m = params.bin.length_m;
params.aeration.nParallelTubes = 12;
params.aeration.bedSideGaugePressure_Pa = 0.0;
params.aeration.bedPorosity = 0.80;
params.aeration.bedParticleDiameter_m = [];
params.aeration.representativeTubeCenterHeight_m = NaN;

params.wire.gauge = '31 AWG';
params.wire.diameter_m = 2.27e-4;
params.wire.rho20_Ohm_m = 1.09e-6;
params.wire.tempCoeff_1_K = 1.7e-4;

params.electrical = struct();
params.electrical.topologyMode = 'paired_series_strings';
params.electrical.manualStringTubeCounts = [];

params.parallel = struct();
params.parallel.enabled = false;
params.parallel.workers = 0;
params.parallel.autoStartPool = true;

params.sweep.voltage_V = 2:2:120;
params.sweep.totalFlow_Lpm = 40:20:400;
params.sweep.pitch_mm = [4 6 8 10 12];
params.sweep.wireDiameter_mm = [0.18 0.20 0.23 0.26 0.29 0.32 0.36 0.41 0.51 0.64];
params.sweep.coilMeanD_mm = [];
params.sweep.wireLength_m = [];
params.sweep.targetDuty = 0.75;
params.sweep.tradeStudyVoltageHalfWindow = 3;
params.sweep.tradeStudyFlowHalfWindow = 2;
params.sweep.tradeStudyFlowPointCount = 100;

params.limits.maxTotalCurrent_A = 15;
params.limits.maxWireTemp_C = 300;
params.limits.maxAirOutletTemp_C = 75;
params.limits.maxPressureDrop_Pa = 1500;
params.limits.maxWaterLoss_kg_day = 5;
params.limits.maxHoleVelocity_m_s = 10;
params.calibration.hWireMultiplier = 1.0;
params.calibration.hWallMultiplier = 1.0;
params.calibration.bedHTMultiplier = 1.0;
params.calibration.evaporationHMultiplier = 1.0;

params.perforation.holesPerTube = 60;
params.perforation.holeDiameter_m = 0.005;
params.perforation.dischargeCoefficient = 0.60;

params.evaporation = struct();
params.evaporation.relativeHumidity = 0.80;
params.evaporation.wettedAreaBasis = 'table3_particle_surface';
params.evaporation.wettedAreaFraction = 1.0;
params.evaporation.accessibleSurfaceFraction = 1.0;
params.evaporation.table3ParticleSurface = defaultTable3ParticleSurface();
params.evaporation.wettedArea_m2 = NaN;
params.evaporation.surfaceTemp_C = params.bin.setpoint_C;
params.evaporation.lewisFactor = 1.0;

params.options.makePlots = true;
params.options.includeTransient = false;
params.options.printSummary = true;

if nargin >= 1 && ~isempty(varargin{1})
    params = mergeStructs(params, varargin{1});
end
params = finalizeDerivedParams(params);

results = runModel(params);
if ~isfield(params, 'options') || ~isfield(params.options, 'printSummary') || ...
        params.options.printSummary
    printSummary(results);
end

if params.options.makePlots
    makePlots(results);
end

if nargout == 0
    assignin('base', 'vermiHeaterResults', results);
end
end

function results = runModel(params)
operatingSweep = runOperatingSweep(params);
results = runScenario(params, params.options.includeTransient, operatingSweep);
results.configurationComparisons = evaluateConfigurationComparisons(params, operatingSweep);
results.activeConfigurationIndex = findActiveConfigurationIndex( ...
    params, results.configurationComparisons);
if isfinite(results.activeConfigurationIndex)
    % Carry the fully evaluated active scenario into the comparison set so
    % all downstream output references the same recommendation source.
    results.configurationComparisons(results.activeConfigurationIndex).scenario = results;
end
results.wireDiameterSweep = evaluateWireDiameterSweep( ...
    params, results.configurationComparisons);
results.coilDiameterSweep = evaluateCoilDiameterSweep( ...
    params, results.configurationComparisons);
results.wireLengthSweep = evaluateWireLengthSweep( ...
    params, results.configurationComparisons);
end

function results = runScenario(params, includeTransient, operatingSweep)
if nargin < 3 || isempty(operatingSweep)
    operatingSweep = runOperatingSweep(params);
end

binModel = buildBinModel(params);
lumpedLoss = lumpedLossModel(params, binModel);
Qreq_W = lumpedLoss.requiredHeat_W;

flowMin100_Lpm = minimumTotalFlowForOutletLimit(params, Qreq_W, 1.0);
flowMinTargetDuty_Lpm = minimumTotalFlowForOutletLimit( ...
    params, Qreq_W, params.sweep.targetDuty);

designPoints = applyScenarioToSweep(operatingSweep, binModel, Qreq_W, params);
recommended = selectRecommendedPoint(designPoints, params);
bestAvailable = selectBestAvailablePoint(designPoints, params);
referencePoint = recommended;
if isempty(fieldnames(referencePoint))
    referencePoint = bestAvailable;
end

recommendedOp = struct();
if ~isempty(fieldnames(referencePoint))
    opParams = params;
    opParams.heaterTube.pitch_m = referencePoint.pitch_mm / 1000;
    recommendedOp = evaluateOperatingPoint( ...
        opParams, referencePoint.voltage_V, referencePoint.totalFlow_Lpm, ...
        params.bin.setpoint_C, params.environment.greenhouseAir_C);

    if includeTransient
        transient = simulateBinPID(opParams, binModel, referencePoint);
    else
        transient = struct();
    end
else
    transient = struct();
end

results = struct();
results.params = params;
results.binModel = binModel;
results.lumpedLoss = lumpedLoss;
results.requiredHeat_W = Qreq_W;
results.minimumFlowAt100Duty_Lpm = flowMin100_Lpm;
results.minimumFlowAtTargetDuty_Lpm = flowMinTargetDuty_Lpm;
results.designPoints = designPoints;
results.recommended = recommended;
results.bestAvailable = bestAvailable;
results.referencePoint = referencePoint;
results.recommendedOperatingPoint = recommendedOp;
results.transient = transient;
end

function comparisons = evaluateConfigurationComparisons(params, operatingSweep)
configs = vesselConfigurationSet(params);
comparisons = struct([]);

for i = 1:numel(configs)
    configParams = applyVesselConfiguration(params, configs(i));
    scenario = runScenario(configParams, false, operatingSweep);

    comparisons(i).id = configs(i).id;
    comparisons(i).label = configs(i).label;
    comparisons(i).color = configs(i).color;
    comparisons(i).scenario = scenario;
end
end

function sweep = evaluateWireDiameterSweep(params, comparisonResults)
sweep = evaluateGeometryTradeSweep( ...
    params, comparisonResults, 'wire_diameter', ...
    params.sweep.wireDiameter_mm(:).', 'diameter_mm', ...
    'Wire diameter', 'mm');
end

function sweep = evaluateCoilDiameterSweep(params, comparisonResults)
values_mm = [];
if isfield(params.sweep, 'coilMeanD_mm')
    values_mm = params.sweep.coilMeanD_mm(:).';
end
sweep = evaluateGeometryTradeSweep( ...
    params, comparisonResults, 'coil_diameter', values_mm, ...
    'coilMeanD_mm', 'Coil mean diameter', 'mm');
end

function sweep = evaluateWireLengthSweep(params, comparisonResults)
values_m = [];
if isfield(params.sweep, 'wireLength_m')
    values_m = params.sweep.wireLength_m(:).';
end
sweep = evaluateGeometryTradeSweep( ...
    params, comparisonResults, 'wire_length', values_m, ...
    'wireLength_m', 'Requested wire length', 'm');
end

function sweep = evaluateGeometryTradeSweep( ...
        params, comparisonResults, studyType, axisValues, axisField, axisLabel, axisUnit)
if nargin < 2 || isempty(comparisonResults)
    operatingSweep = runOperatingSweep(params);
    comparisonResults = evaluateConfigurationComparisons(params, operatingSweep);
end

if isempty(axisValues)
    sweep = struct([]);
    return;
end

configs = vesselConfigurationSet(params);
sweep = struct([]);
selectedCriteriaKeys = cell(numel(configs), 1);
availableCriteriaKeys = cell(numel(configs), 1);
feasibleCriteriaKeys = cell(numel(configs), 1);

for iCfg = 1:numel(configs)
    targetRec = struct();
    if iCfg <= numel(comparisonResults) && ...
            isfield(comparisonResults(iCfg), 'scenario')
        targetRec = scenarioReferencePoint(comparisonResults(iCfg).scenario);
    end

    sweep(iCfg).studyType = studyType;
    sweep(iCfg).xField = axisField;
    sweep(iCfg).xLabel = axisLabel;
    sweep(iCfg).xUnit = axisUnit;
    sweep(iCfg).xValues = axisValues;
    sweep(iCfg).id = configs(iCfg).id;
    sweep(iCfg).label = configs(iCfg).label;
    sweep(iCfg).color = configs(iCfg).color;
    sweep(iCfg).(axisField) = axisValues;
    sweep(iCfg).targetPower_W = NaN;
    sweep(iCfg).targetFlow_Lpm = NaN;
    sweep(iCfg).targetCurrent_A = NaN;
    sweep(iCfg).targetVoltage_V = NaN;
    sweep(iCfg).targetPitch_mm = NaN;
    sweep(iCfg).optimalPower_W = nan(size(axisValues));
    sweep(iCfg).optimalFlow_Lpm = nan(size(axisValues));
    sweep(iCfg).optimalCurrent_A = nan(size(axisValues));
    sweep(iCfg).optimalVoltage_V = nan(size(axisValues));
    sweep(iCfg).optimalPitch_mm = nan(size(axisValues));
    sweep(iCfg).optimalDuty = nan(size(axisValues));
    sweep(iCfg).optimalQtoBed_W = nan(size(axisValues));
    sweep(iCfg).optimalWireMax_C = nan(size(axisValues));
    sweep(iCfg).optimalColderNode_C = nan(size(axisValues));
    sweep(iCfg).selectionScore = nan(size(axisValues));
    sweep(iCfg).selectedPower_W = nan(size(axisValues));
    sweep(iCfg).selectedFlow_Lpm = nan(size(axisValues));
    sweep(iCfg).selectedCurrent_A = nan(size(axisValues));
    sweep(iCfg).selectedVoltage_V = nan(size(axisValues));
    sweep(iCfg).selectedPitch_mm = nan(size(axisValues));
    sweep(iCfg).selectedDuty = nan(size(axisValues));
    sweep(iCfg).selectedQtoBed_W = nan(size(axisValues));
    sweep(iCfg).selectedWireMax_C = nan(size(axisValues));
    sweep(iCfg).selectedColderNode_C = nan(size(axisValues));
    sweep(iCfg).selectedSelectionScore = nan(size(axisValues));
    sweep(iCfg).selectionMode = repmat({''}, size(axisValues));
    sweep(iCfg).isFeasible = false(size(axisValues));
    sweep(iCfg).isConstraintSafe = false(size(axisValues));
    sweep(iCfg).bestAvailablePower_W = nan(size(axisValues));
    sweep(iCfg).bestAvailableFlow_Lpm = nan(size(axisValues));
    sweep(iCfg).bestAvailableCurrent_A = nan(size(axisValues));
    sweep(iCfg).bestAvailableVoltage_V = nan(size(axisValues));
    sweep(iCfg).bestAvailablePitch_mm = nan(size(axisValues));
    sweep(iCfg).bestAvailableDuty = nan(size(axisValues));
    sweep(iCfg).bestAvailableQtoBed_W = nan(size(axisValues));
    sweep(iCfg).bestAvailableWireMax_C = nan(size(axisValues));
    sweep(iCfg).bestAvailableColderNode_C = nan(size(axisValues));
    sweep(iCfg).bestAvailableMode = repmat({''}, size(axisValues));
    sweep(iCfg).bestIndex = NaN;
    sweep(iCfg).bestAvailableIndex = NaN;
    selectedCriteriaKeys{iCfg} = inf(numel(axisValues), 14);
    availableCriteriaKeys{iCfg} = inf(numel(axisValues), 14);
    feasibleCriteriaKeys{iCfg} = inf(numel(axisValues), 14);

    if ~isempty(fieldnames(targetRec))
        sweep(iCfg).targetPower_W = targetRec.totalPower_W;
        sweep(iCfg).targetFlow_Lpm = targetRec.totalFlow_Lpm;
        sweep(iCfg).targetCurrent_A = targetRec.totalCurrent_A;
        sweep(iCfg).targetVoltage_V = targetRec.voltage_V;
        sweep(iCfg).targetPitch_mm = targetRec.pitch_mm;
    end
end

for iVal = 1:numel(axisValues)
    trialParams = applyGeometryTradeValue(params, studyType, axisValues(iVal));

    for iCfg = 1:numel(configs)
        targetRec = struct();
        if isfinite(sweep(iCfg).targetPower_W)
            targetRec.totalPower_W = sweep(iCfg).targetPower_W;
            targetRec.totalFlow_Lpm = sweep(iCfg).targetFlow_Lpm;
            targetRec.totalCurrent_A = sweep(iCfg).targetCurrent_A;
            targetRec.voltage_V = sweep(iCfg).targetVoltage_V;
            targetRec.pitch_mm = sweep(iCfg).targetPitch_mm;
        end

        configTrialParams = applyVesselConfiguration(trialParams, configs(iCfg));
        operatingSweep = runTargetedOperatingSweep(configTrialParams, targetRec);
        binModel = buildBinModel(configTrialParams);
        lumpedLoss = lumpedLossModel(configTrialParams, binModel);
        designPoints = applyScenarioToSweep( ...
            operatingSweep, binModel, lumpedLoss.requiredHeat_W, configTrialParams);
        feasible = [designPoints.isFeasible];
        constraintSafe = [designPoints.isConstraintSafe];
        [bestAvailable, bestAvailableCriteria] = selectBestPointByCriteria( ...
            designPoints, configTrialParams);
        bestAvailableMode = criteriaSelectionMode(bestAvailable, configTrialParams);

        if ~isempty(fieldnames(bestAvailable))
            sweep(iCfg).bestAvailablePower_W(iVal) = bestAvailable.totalPower_W;
            sweep(iCfg).bestAvailableFlow_Lpm(iVal) = bestAvailable.totalFlow_Lpm;
            sweep(iCfg).bestAvailableCurrent_A(iVal) = bestAvailable.totalCurrent_A;
            sweep(iCfg).bestAvailableVoltage_V(iVal) = bestAvailable.voltage_V;
            sweep(iCfg).bestAvailablePitch_mm(iVal) = bestAvailable.pitch_mm;
            sweep(iCfg).bestAvailableDuty(iVal) = bestAvailable.dutyNeeded;
            sweep(iCfg).bestAvailableQtoBed_W(iVal) = bestAvailable.QtoBed_W;
            sweep(iCfg).bestAvailableWireMax_C(iVal) = bestAvailable.wireMax_C;
            sweep(iCfg).bestAvailableColderNode_C(iVal) = min( ...
                bestAvailable.TbottomFullPower_C, bestAvailable.TtopFullPower_C);
            sweep(iCfg).bestAvailableMode{iVal} = bestAvailableMode;
            availableCriteriaKeys{iCfg}(iVal, :) = bestAvailableCriteria;
        end

        if any(constraintSafe)
            [best, bestCriteria] = selectBestPointByCriteria( ...
                designPoints(constraintSafe), configTrialParams);
            selectionMode = criteriaSelectionMode(best, configTrialParams);
            sweep(iCfg).selectedPower_W(iVal) = best.totalPower_W;
            sweep(iCfg).selectedFlow_Lpm(iVal) = best.totalFlow_Lpm;
            sweep(iCfg).selectedCurrent_A(iVal) = best.totalCurrent_A;
            sweep(iCfg).selectedVoltage_V(iVal) = best.voltage_V;
            sweep(iCfg).selectedPitch_mm(iVal) = best.pitch_mm;
            sweep(iCfg).selectedDuty(iVal) = best.dutyNeeded;
            sweep(iCfg).selectedQtoBed_W(iVal) = best.QtoBed_W;
            sweep(iCfg).selectedWireMax_C(iVal) = best.wireMax_C;
            sweep(iCfg).selectedColderNode_C(iVal) = min( ...
                best.TbottomFullPower_C, best.TtopFullPower_C);
            sweep(iCfg).selectedSelectionScore(iVal) = selectionStageCode(selectionMode);
            sweep(iCfg).selectionMode{iVal} = selectionMode;
            sweep(iCfg).isConstraintSafe(iVal) = true;
            selectedCriteriaKeys{iCfg}(iVal, :) = bestCriteria;
        end

        if any(feasible)
            [bestFeasible, feasibleCriteria] = selectBestPointByCriteria( ...
                designPoints(feasible), configTrialParams);
            sweep(iCfg).optimalPower_W(iVal) = bestFeasible.totalPower_W;
            sweep(iCfg).optimalFlow_Lpm(iVal) = bestFeasible.totalFlow_Lpm;
            sweep(iCfg).optimalCurrent_A(iVal) = bestFeasible.totalCurrent_A;
            sweep(iCfg).optimalVoltage_V(iVal) = bestFeasible.voltage_V;
            sweep(iCfg).optimalPitch_mm(iVal) = bestFeasible.pitch_mm;
            sweep(iCfg).optimalDuty(iVal) = bestFeasible.dutyNeeded;
            sweep(iCfg).optimalQtoBed_W(iVal) = bestFeasible.QtoBed_W;
            sweep(iCfg).optimalWireMax_C(iVal) = bestFeasible.wireMax_C;
            sweep(iCfg).optimalColderNode_C(iVal) = min( ...
                bestFeasible.TbottomFullPower_C, bestFeasible.TtopFullPower_C);
            sweep(iCfg).selectionScore(iVal) = 1;
            sweep(iCfg).isFeasible(iVal) = true;
            feasibleCriteriaKeys{iCfg}(iVal, :) = feasibleCriteria;
        end
    end
end

for iCfg = 1:numel(configs)
    validIdx = find(sweep(iCfg).isConstraintSafe);
    if ~isempty(validIdx)
        [~, order] = sortrows(selectedCriteriaKeys{iCfg}(validIdx, :));
        sweep(iCfg).bestIndex = validIdx(order(1));
    end
    validAvailableIdx = find(all(isfinite(availableCriteriaKeys{iCfg}), 2));
    if ~isempty(validAvailableIdx)
        [~, order] = sortrows(availableCriteriaKeys{iCfg}(validAvailableIdx, :));
        sweep(iCfg).bestAvailableIndex = validAvailableIdx(order(1));
    end
end
end

function trialParams = applyGeometryTradeValue(params, studyType, value)
trialParams = params;
trialParams.heaterTube.manualWireLength_m = NaN;
switch studyType
    case 'wire_diameter'
        trialParams.wire.diameter_m = value / 1000;
        trialParams.wire.gauge = nearestAwgLabel(trialParams.wire.diameter_m);
    case 'coil_diameter'
        trialParams.heaterTube.coilMeanD_m = value / 1000;
    case 'wire_length'
        trialParams.heaterTube.manualWireLength_m = value;
end
end

function configs = vesselConfigurationSet(params)
ventDia_m = params.bin.ventHoleDiameter_m;
if ventDia_m <= 0
    ventDia_m = params.bin.comparisonVentHoleDiameter_m;
end

configs = struct([]);
configs(1).id = 'covered';
configs(1).label = 'Covered top';
configs(1).openTop = false;
configs(1).ventHoleCount = 0;
configs(1).ventHoleDiameter_m = 0.0;
configs(1).color = [0.00 0.45 0.74];

configs(2).id = 'covered_vent';
configs(2).label = 'Covered top + vent';
configs(2).openTop = false;
configs(2).ventHoleCount = max(1, params.bin.ventHoleCount);
configs(2).ventHoleDiameter_m = ventDia_m;
configs(2).color = [0.85 0.33 0.10];

configs(3).id = 'open';
configs(3).label = 'Uncovered top';
configs(3).openTop = true;
configs(3).ventHoleCount = 0;
configs(3).ventHoleDiameter_m = 0.0;
configs(3).color = [0.47 0.67 0.19];
end

function cfgParams = applyVesselConfiguration(params, config)
cfgParams = params;
cfgParams.bin.openTop = config.openTop;
cfgParams.bin.ventHoleCount = config.ventHoleCount;
cfgParams.bin.ventHoleDiameter_m = config.ventHoleDiameter_m;
end

function idx = findActiveConfigurationIndex(params, comparisons)
idx = NaN;
for i = 1:numel(comparisons)
    cfgParams = comparisons(i).scenario.params;
    sameOpenTop = logical(cfgParams.bin.openTop) == logical(params.bin.openTop);
    sameVentCount = cfgParams.bin.ventHoleCount == params.bin.ventHoleCount;
    sameVentDia = abs(cfgParams.bin.ventHoleDiameter_m - ...
        params.bin.ventHoleDiameter_m) < 1e-12;
    if sameOpenTop && sameVentCount && sameVentDia
        idx = i;
        return;
    end
end
end

function [activeCfg, activeScenario, rec] = comparisonBasedReference(results)
activeCfg = struct();
activeScenario = results;
rec = scenarioReferencePoint(results);

if isfield(results, 'activeConfigurationIndex') && ...
        isfinite(results.activeConfigurationIndex) && ...
        results.activeConfigurationIndex >= 1 && ...
        results.activeConfigurationIndex <= numel(results.configurationComparisons)
    activeCfg = results.configurationComparisons(results.activeConfigurationIndex);
    activeScenario = activeCfg.scenario;
    rec = scenarioReferencePoint(activeScenario);
end
end

function rec = scenarioReferencePoint(scenario)
rec = struct();
if isfield(scenario, 'recommended') && ~isempty(fieldnames(scenario.recommended))
    rec = scenario.recommended;
elseif isfield(scenario, 'referencePoint') && ~isempty(fieldnames(scenario.referencePoint))
    rec = scenario.referencePoint;
elseif isfield(scenario, 'bestAvailable') && ~isempty(fieldnames(scenario.bestAvailable))
    rec = scenario.bestAvailable;
end
end

function [best, bestCriteria] = selectBestPointByCriteria(candidates, params)
if isempty(candidates)
    best = struct();
    bestCriteria = inf(1, 14);
    return;
end

criteria = criteriaMatrix(candidates, params);
[~, order] = sortrows(criteria);
best = candidates(order(1));
bestCriteria = criteria(order(1), :);
end

function criteria = criteriaMatrix(points, params)
criteria = inf(numel(points), 14);
for i = 1:numel(points)
    criteria(i, :) = pointCriteriaRow(points(i), params);
end
end

function row = pointCriteriaRow(pt, params)
Tbottom_C = pt.TbottomFullPower_C;
Ttop_C = pt.TtopFullPower_C;
if ~isfinite(Tbottom_C)
    Tbottom_C = -1e6;
end
if ~isfinite(Ttop_C)
    Ttop_C = -1e6;
end
spread_C = abs(Tbottom_C - Ttop_C);

% Lexicographic order: hard constraints first, then maximize the colder
% bed-node temperature, then minimize power, then minimize spread.
row = [ ...
    max(0, pt.totalCurrent_A - params.limits.maxTotalCurrent_A), ...
    max(0, pt.wireMax_C - params.limits.maxWireTemp_C), ...
    max(0, pt.airOutlet_C - params.limits.maxAirOutletTemp_C), ...
    max(0, pt.deltaP_Pa - params.limits.maxPressureDrop_Pa), ...
    max(0, pt.waterLoss_kg_day - params.limits.maxWaterLoss_kg_day), ...
    max(0, pt.holeVelocity_m_s - params.limits.maxHoleVelocity_m_s), ...
    max(0, params.optimization.minBottomTemp_C - Tbottom_C), ...
    max(0, params.optimization.minTopTemp_C - Ttop_C), ...
    max(0, spread_C - params.optimization.spreadLimit_C), ...
    -min(Tbottom_C, Ttop_C), ...
    pt.totalPower_W, ...
    spread_C, ...
    -0.5 * (Tbottom_C + Ttop_C), ...
    max(0, pt.dutyNeeded - 1.0)];
end

function mode = criteriaSelectionMode(pt, params)
if isempty(fieldnames(pt))
    mode = '';
elseif pt.isFeasible
    mode = 'feasible';
elseif isfield(pt, 'isConstraintSafe') && pt.isConstraintSafe
    mode = 'fallback';
else
    mode = 'last_resort';
end
end

function code = selectionStageCode(mode)
switch mode
    case 'feasible'
        code = 1;
    case 'fallback'
        code = 2;
    otherwise
        code = 3;
end
end

function params = defaultParams()
params = struct();

params.environment = struct();
params.environment.greenhouseAir_C = -10;
params.environment.pressure_Pa = 101325;
params.environment.referencePlateLength_m = 0.61;

params.bin = struct();
params.bin.length_m = 1.22;
params.bin.width_m = 0.61;
params.bin.height_m = 0.61;
params.bin.openTop = true;
params.bin.fillFraction = 0.80;
params.bin.bulkDensity_kg_m3 = 650;
params.bin.cp_J_kgK = 3200;
params.bin.k_W_mK = 0.45;
params.bin.ventHoleCount = 1;
params.bin.ventHoleDiameter_m = 0.0;
params.bin.comparisonVentHoleDiameter_m = 0.025;
params.bin.setpoint_C = 20;
params.bin.minSafe_C = 15;
params.bin.maxSafe_C = 25;
params.bin.heatToBottomFraction = 0.60;

params.binWall = struct();
params.binWall.sheetThickness_m = 1.5e-3;
params.binWall.sheetK_W_mK = 50;
params.binWall.insulationThickness_m = 0.0508;
params.binWall.insulationK_W_mK = 0.029;
params.binWall.internalH_W_m2K = 6;

params.heaterTube = struct();
params.heaterTube.length_m = 0.30;
params.heaterTube.ID_m = 0.015;
params.heaterTube.OD_m = 0.019;
params.heaterTube.wallK_W_mK = 16;
params.heaterTube.insulationThickness_m = 0.010;
params.heaterTube.insulationK_W_mK = 0.040;
params.heaterTube.externalH_W_m2K = 10;
params.heaterTube.coilMeanD_m = 0.0115;
params.heaterTube.pitch_m = 0.008;
params.heaterTube.coilSpan_m = 0.30;
params.heaterTube.manualWireLength_m = NaN;
params.heaterTube.segments = 50;
params.heaterTube.extraMinorK = 1.0;

params.aeration = struct();
params.aeration.length_m = 1.22;
params.aeration.ID_m = 0.040;
params.aeration.OD_m = 0.044;
params.aeration.wallK_W_mK = 16;
params.aeration.segments = 40;
params.aeration.bedH_W_m2K = 8;
params.aeration.nParallelTubes = 12;
params.aeration.headerEnabled = true;
params.aeration.headerID_m = 0.050;
params.aeration.headerLossMultiplier = 1.0;
params.aeration.splitterOutletCount = 4;
params.aeration.splitterInletID_m = 0.019;
params.aeration.branchConnectorID_m = 0.019;
params.aeration.branchConnectorLength_m = 0.0;
params.aeration.splitterBodyLossMultiplier = 1.0;
params.aeration.contractionModel = 'idelchik_conical_bellmouth_without_end_wall';
params.aeration.contractionConeAngle_deg = 60.0;
params.aeration.bedSideGaugePressure_Pa = 0.0;
params.aeration.bedPorosity = 0.80;
params.aeration.bedParticleDiameter_m = [];
params.aeration.representativeTubeCenterHeight_m = NaN;

params.wire = struct();
params.wire.name = 'Nichrome 80/20';
params.wire.gauge = '31 AWG';
params.wire.diameter_m = 2.27e-4;
params.wire.rho20_Ohm_m = 1.09e-6;
params.wire.tempCoeff_1_K = 1.7e-4;

params.sweep = struct();
params.sweep.voltage_V = 2:2:120;
params.sweep.totalFlow_Lpm = 40:20:400;
params.sweep.pitch_mm = [4 6 8 10 12];
params.sweep.wireDiameter_mm = [0.18 0.20 0.23 0.26 0.29 0.32 0.36 0.41 0.51 0.64];
params.sweep.coilMeanD_mm = [];
params.sweep.wireLength_m = [];
params.sweep.targetDuty = 0.75;
params.sweep.tradeStudyVoltageHalfWindow = 3;
params.sweep.tradeStudyFlowHalfWindow = 2;
params.sweep.tradeStudyFlowPointCount = 100;

params.limits = struct();
params.limits.maxTotalCurrent_A = 15;
params.limits.maxWireTemp_C = 300;
params.limits.maxAirOutletTemp_C = 75;
params.limits.maxPressureDrop_Pa = 1500;
params.limits.maxWaterLoss_kg_day = 5;
params.limits.maxHoleVelocity_m_s = 10;

params.optimization = struct();
params.optimization.spreadLimit_C = 10;
params.optimization.minBottomTemp_C = params.bin.minSafe_C;
params.optimization.minTopTemp_C = params.bin.minSafe_C;

params.control = struct();
params.control.duration_h = 72;
params.control.dt_s = 60;
params.control.initialBottom_C = 12;
params.control.initialTop_C = 12;
params.control.Kp = 4.5;
params.control.Ki = 0.003;
params.control.Kd = 0.0;
params.control.highCutout_C = 28;

params.perforation = struct();
params.perforation.holesPerTube = 60;
params.perforation.holeDiameter_m = 0.005;
params.perforation.dischargeCoefficient = 0.60;

params.evaporation = struct();
params.evaporation.relativeHumidity = 0.80;
params.evaporation.wettedAreaBasis = 'table3_particle_surface';
params.evaporation.wettedAreaFraction = 1.0;
params.evaporation.accessibleSurfaceFraction = 1.0;
params.evaporation.table3ParticleSurface = defaultTable3ParticleSurface();
params.evaporation.wettedArea_m2 = NaN;
params.evaporation.surfaceTemp_C = params.bin.setpoint_C;
params.evaporation.lewisFactor = 1.0;

params.calibration = struct();
params.calibration.hWireMultiplier = 1.0;
params.calibration.hWallMultiplier = 1.0;
params.calibration.bedHTMultiplier = 1.0;
params.calibration.evaporationHMultiplier = 1.0;

params.options = struct();
params.options.makePlots = true;
params.options.includeTransient = false;
params.options.printSummary = true;
params = finalizeDerivedParams(params);
end

function dst = mergeStructs(dst, src)
if isempty(src)
    return;
end

fields = fieldnames(src);
for i = 1:numel(fields)
    name = fields{i};
    srcVal = src.(name);
    if isstruct(srcVal) && isfield(dst, name) && isstruct(dst.(name))
        dst.(name) = mergeStructs(dst.(name), srcVal);
    else
        dst.(name) = srcVal;
    end
end
end

function params = finalizeDerivedParams(params)
if isfield(params, 'aeration') && isstruct(params.aeration) && ...
        isfield(params.aeration, 'nParallelTubes') && ...
        ~isempty(params.aeration.nParallelTubes)
    params.aeration.splitterOutletCount = max(1, round(params.aeration.nParallelTubes));
end

if ~isfield(params, 'electrical') || ~isstruct(params.electrical)
    params.electrical = struct();
end
if ~isfield(params.electrical, 'topologyMode') || ...
        isempty(params.electrical.topologyMode)
    params.electrical.topologyMode = 'paired_series_strings';
end
if ~isfield(params.electrical, 'manualStringTubeCounts') || ...
        isempty(params.electrical.manualStringTubeCounts)
    params.electrical.manualStringTubeCounts = [];
end
params.electrical.stringTubeCounts = preferredElectricalStringTubeCounts( ...
    max(1, round(params.aeration.nParallelTubes)), params.electrical);
params.electrical.topologyLabel = electricalTopologyLabel( ...
    params.electrical.stringTubeCounts);

if ~isfield(params, 'parallel') || ~isstruct(params.parallel)
    params.parallel = struct();
end
if ~isfield(params.parallel, 'enabled') || isempty(params.parallel.enabled)
    params.parallel.enabled = false;
end
if ~isfield(params.parallel, 'workers') || isempty(params.parallel.workers)
    params.parallel.workers = 0;
end
if ~isfield(params.parallel, 'autoStartPool') || isempty(params.parallel.autoStartPool)
    params.parallel.autoStartPool = true;
end

if ~isfield(params, 'evaporation') || ~isstruct(params.evaporation)
    return;
end

basis = 'top_plan_area';
if isfield(params.evaporation, 'wettedAreaBasis') && ...
        ~isempty(params.evaporation.wettedAreaBasis)
    basis = lower(string(params.evaporation.wettedAreaBasis));
end

wetFraction = 1.0;
if isfield(params.evaporation, 'wettedAreaFraction') && ...
        ~isempty(params.evaporation.wettedAreaFraction)
    wetFraction = max(0.0, params.evaporation.wettedAreaFraction);
end

accessibleFraction = wetFraction;
if isfield(params.evaporation, 'accessibleSurfaceFraction') && ...
        ~isempty(params.evaporation.accessibleSurfaceFraction)
    accessibleFraction = max(0.0, params.evaporation.accessibleSurfaceFraction);
end

switch basis
    case "manual"
        if ~isfield(params.evaporation, 'wettedArea_m2') || ...
                isempty(params.evaporation.wettedArea_m2)
            params.evaporation.wettedArea_m2 = derivedWettedArea_m2(params, wetFraction);
        end
    case "table3_particle_surface"
        params.evaporation.wettedArea_m2 = derivedWettedAreaParticleSurface_m2( ...
            params, accessibleFraction);
        params.evaporation.wettedAreaBasis = char(basis);
        params.evaporation.accessibleSurfaceFraction = accessibleFraction;
    otherwise
        params.evaporation.wettedArea_m2 = derivedWettedArea_m2(params, wetFraction);
        params.evaporation.wettedAreaBasis = char(basis);
        params.evaporation.wettedAreaFraction = wetFraction;
end
end

function counts = preferredElectricalStringTubeCounts(nTubes, electrical)
nTubes = max(1, round(nTubes));
counts = [];

if isfield(electrical, 'manualStringTubeCounts') && ...
        ~isempty(electrical.manualStringTubeCounts)
    manual = reshape(electrical.manualStringTubeCounts, 1, []);
    manual = round(manual(isfinite(manual) & manual > 0));
    if ~isempty(manual) && sum(manual) == nTubes
        counts = manual;
        return;
    end
end

mode = 'paired_series_strings';
if isfield(electrical, 'topologyMode') && ~isempty(electrical.topologyMode)
    mode = lower(string(electrical.topologyMode));
end

switch mode
    case "all_parallel"
        counts = ones(1, nTubes);
    otherwise
        if nTubes == 1
            counts = 1;
        elseif mod(nTubes, 2) == 0
            counts = 2 * ones(1, nTubes / 2);
        else
            % Odd tube counts are handled by adding the extra tube in
            % series with one of the otherwise paired strings.
            counts = [3, 2 * ones(1, (nTubes - 3) / 2)];
        end
end
end

function label = electricalTopologyLabel(counts)
counts = reshape(round(counts(isfinite(counts) & counts > 0)), 1, []);
if isempty(counts)
    label = 'unspecified';
    return;
end

uniqueCounts = unique(counts);
parts = strings(1, numel(uniqueCounts));
for i = 1:numel(uniqueCounts)
    nSeries = uniqueCounts(i);
    nStrings = sum(counts == nSeries);
    parts(i) = sprintf('%dx%dS', nStrings, nSeries);
end
label = strjoin(cellstr(parts), ' + ');
end

function area_m2 = derivedWettedArea_m2(params, wetFraction)
% [R5] Macroscopic wetted area approximation for the evaporation model:
% A_wet = f_wet * L * W
area_m2 = max(0.0, wetFraction) * params.bin.length_m * params.bin.width_m;
end

function table3 = defaultTable3ParticleSurface()
table3 = struct();
table3.blendMode = 'mean_of_windrow_and_vermicompost';
table3.representativeDiameters_mm = struct( ...
    'gt40', 60.0, ...
    'mm20to40', 30.0, ...
    'mm10to20', 15.0, ...
    'lt10', 5.0);
table3.windrowFractions_pct = struct( ...
    'gt40', 9.8, ...
    'mm20to40', 26.6, ...
    'mm10to20', 26.6, ...
    'lt10', 37.0);
table3.vermicompostFractions_pct = struct( ...
    'gt40', 7.9, ...
    'mm20to40', 6.9, ...
    'mm10to20', 19.8, ...
    'lt10', 65.3);
end

function area_m2 = derivedWettedAreaParticleSurface_m2(params, accessibleFraction)
% [R5,R11] Particle-surface wetted-area estimate from the Table 3 particle
% size fractions. Each size class is represented by an equivalent-sphere
% diameter d_i, the blended mass fraction is w_i, the Sauter mean diameter
% is d32 = 1/sum(w_i/d_i), and the corresponding external area is
% A_wet = f_access * 6 * V_fill / d32.
d32_m = table3ParticleSauterDiameter_m(params);
fillVolume_m3 = params.bin.length_m * params.bin.width_m * ...
    params.bin.height_m * params.bin.fillFraction;
area_m2 = max(0.0, accessibleFraction) * 6 * fillVolume_m3 / max(d32_m, 1e-12);
end

function d32_m = table3ParticleSauterDiameter_m(params)
d32_m = NaN;
table3 = defaultTable3ParticleSurface();
if isfield(params, 'evaporation') && isstruct(params.evaporation) && ...
        isfield(params.evaporation, 'table3ParticleSurface') && ...
        isstruct(params.evaporation.table3ParticleSurface)
    table3 = mergeStructs(table3, params.evaporation.table3ParticleSurface);
end

[weights, diameters_m] = table3ParticleBlend(table3);
if isempty(weights) || any(~isfinite(diameters_m))
    return;
end

invD32_m_inv = sum(weights ./ max(diameters_m, 1e-12));
d32_m = 1 / max(invD32_m_inv, 1e-12);
end

function [weights, diameters_m] = table3ParticleBlend(table3)
keys = {'gt40', 'mm20to40', 'mm10to20', 'lt10'};
weights = zeros(1, numel(keys));
diameters_m = zeros(1, numel(keys));

for i = 1:numel(keys)
    key = keys{i};
    diameters_m(i) = table3.representativeDiameters_mm.(key) / 1000;
end

blendMode = 'mean_of_windrow_and_vermicompost';
if isfield(table3, 'blendMode') && ~isempty(table3.blendMode)
    blendMode = lower(string(table3.blendMode));
end

switch blendMode
    case "mean_of_windrow_and_vermicompost"
        for i = 1:numel(keys)
            key = keys{i};
            weights(i) = 0.5 * ( ...
                table3.windrowFractions_pct.(key) + ...
                table3.vermicompostFractions_pct.(key));
        end
    otherwise
        for i = 1:numel(keys)
            key = keys{i};
            weights(i) = table3.windrowFractions_pct.(key);
        end
end

weights = max(weights, 0);
weights = weights / max(sum(weights), 1e-12);
end

function binModel = buildBinModel(params)
L = params.bin.length_m;
W = params.bin.width_m;
H = params.bin.height_m;

sideArea_m2 = 2 * H * L + 2 * H * W;
bottomArea_m2 = L * W;
topArea_m2 = L * W;
totalArea_m2 = sideArea_m2 + bottomArea_m2 + topArea_m2;
ventHoleArea_m2 = min(topArea_m2, ...
    max(0, params.bin.ventHoleCount) * pi * ...
    max(params.bin.ventHoleDiameter_m, 0)^2 / 4);
coveredTopArea_m2 = max(topArea_m2 - ventHoleArea_m2, 0);

% [R4] Free convection from a vertical surface in stagnant air.
hExt_W_m2K = naturalConvectionVerticalPlateH( ...
    params.bin.setpoint_C, ...
    params.environment.greenhouseAir_C, ...
    params.environment.referencePlateLength_m, ...
    params.environment.pressure_Pa);

% [R5] Thermal resistances in series: U = 1/sum(R_i)
Uwall_W_m2K = overallWallU(params.binWall, hExt_W_m2K);

AbottomNode_m2 = bottomArea_m2 + 0.5 * sideArea_m2;
topWallArea_m2 = 0.5 * sideArea_m2;

if params.bin.openTop
    % Model the top as a direct bed-to-ambient horizontal free-convection surface.
    topCharLength_m = topArea_m2 / max(2 * (L + W), 1e-9);
    hTopAmbient_W_m2K = naturalConvectionHorizontalPlateH( ...
        params.bin.setpoint_C, ...
        params.environment.greenhouseAir_C, ...
        topCharLength_m, ...
        params.environment.pressure_Pa);
    UAtopAmbient_W_K = Uwall_W_m2K * topWallArea_m2 + ...
        hTopAmbient_W_m2K * topArea_m2;
    UAtotalAmbient_W_K = Uwall_W_m2K * (sideArea_m2 + bottomArea_m2) + ...
        hTopAmbient_W_m2K * topArea_m2;
    topMode = 'open top';
    effectiveOpenTopArea_m2 = topArea_m2;
else
    if ventHoleArea_m2 > 0
        % Extension of the covered/open-top PDF model: treat the vent hole
        % as a locally open area and the remaining top area as covered.
        ventPerimeter_m = max(pi * max(params.bin.ventHoleDiameter_m, 0) * ...
            max(params.bin.ventHoleCount, 0), 1e-9);
        ventCharLength_m = ventHoleArea_m2 / ventPerimeter_m;
        hTopAmbient_W_m2K = naturalConvectionHorizontalPlateH( ...
            params.bin.setpoint_C, ...
            params.environment.greenhouseAir_C, ...
            ventCharLength_m, ...
            params.environment.pressure_Pa);
        UAtopAmbient_W_K = Uwall_W_m2K * (coveredTopArea_m2 + topWallArea_m2) + ...
            hTopAmbient_W_m2K * ventHoleArea_m2;
        UAtotalAmbient_W_K = Uwall_W_m2K * ...
            (sideArea_m2 + bottomArea_m2 + coveredTopArea_m2) + ...
            hTopAmbient_W_m2K * ventHoleArea_m2;
        topMode = 'covered top with vent hole';
        effectiveOpenTopArea_m2 = ventHoleArea_m2;
    else
        hTopAmbient_W_m2K = Uwall_W_m2K;
        UAtopAmbient_W_K = Uwall_W_m2K * (topArea_m2 + topWallArea_m2);
        UAtotalAmbient_W_K = Uwall_W_m2K * totalArea_m2;
        topMode = 'covered top';
        effectiveOpenTopArea_m2 = 0;
    end
end

volume_m3 = L * W * H * params.bin.fillFraction;
mass_kg = volume_m3 * params.bin.bulkDensity_kg_m3;
capacityTotal_J_K = mass_kg * params.bin.cp_J_kgK;

% [R5] Simple internal conductive coupling between bottom and top nodes.
interfaceArea_m2 = L * W;
interfaceThickness_m = max(H / 2, 1e-6);
UAinternal_W_K = params.bin.k_W_mK * interfaceArea_m2 / interfaceThickness_m;

% [R5] Bi = h*Lc/k, with Lc = V/A
Lc_m = volume_m3 / max(totalArea_m2, 1e-9);
Bi_lumped = hExt_W_m2K * Lc_m / max(params.bin.k_W_mK, 1e-9);

binModel = struct();
binModel.totalArea_m2 = totalArea_m2;
binModel.topMode = topMode;
binModel.topArea_m2 = topArea_m2;
binModel.coveredTopArea_m2 = coveredTopArea_m2;
binModel.openTopArea_m2 = effectiveOpenTopArea_m2;
binModel.ventHoleArea_m2 = ventHoleArea_m2;
binModel.externalH_W_m2K = hExt_W_m2K;
binModel.topAmbientH_W_m2K = hTopAmbient_W_m2K;
binModel.Uwall_W_m2K = Uwall_W_m2K;
binModel.UAtotalAmbient_W_K = UAtotalAmbient_W_K;
binModel.UAbottomAmbient_W_K = Uwall_W_m2K * AbottomNode_m2;
binModel.UAtopAmbient_W_K = UAtopAmbient_W_K;
binModel.UAinternal_W_K = UAinternal_W_K;
binModel.Ctotal_J_K = capacityTotal_J_K;
binModel.Cbottom_J_K = 0.5 * capacityTotal_J_K;
binModel.Ctop_J_K = 0.5 * capacityTotal_J_K;
binModel.mass_kg = mass_kg;
binModel.volume_m3 = volume_m3;
binModel.characteristicLength_m = Lc_m;
binModel.Bi_lumped = Bi_lumped;
end

function loss = lumpedLossModel(params, binModel)
Tinf_C = params.environment.greenhouseAir_C;
Tset_C = params.bin.setpoint_C;
UA_W_K = binModel.UAtotalAmbient_W_K;
C_J_K = binModel.Ctotal_J_K;

% [R5] Q_loss = UA * (T_bin - T_inf)
Qreq_W = UA_W_K * (Tset_C - Tinf_C);

% [R5] tau = C / UA
tau_s = C_J_K / max(UA_W_K, 1e-9);

% [R5] T(t) = T_inf + (T0 - T_inf) * exp(-t/tau)
t_h = linspace(0, 72, 145).';
cooldown_C = Tinf_C + (Tset_C - Tinf_C) * exp(-(t_h * 3600) / max(tau_s, 1e-9));

loss = struct();
loss.requiredHeat_W = Qreq_W;
loss.UA_W_K = UA_W_K;
loss.tau_h = tau_s / 3600;
loss.cooldown_time_h = t_h;
loss.cooldown_C = cooldown_C;
loss.Bi_lumped = binModel.Bi_lumped;
loss.lumpedStrictlyValid = binModel.Bi_lumped < 0.1;
end

function U_W_m2K = overallWallU(wall, hOutside_W_m2K)
% [R5] One-dimensional thermal resistance network through the wall.
Rin = 1 / max(wall.internalH_W_m2K, 1e-9);
Rsheet = wall.sheetThickness_m / max(wall.sheetK_W_mK, 1e-9);
Rins = wall.insulationThickness_m / max(wall.insulationK_W_mK, 1e-9);
Rout = 1 / max(hOutside_W_m2K, 1e-9);
U_W_m2K = 1 / (Rin + Rsheet + Rins + Rout);
end

function h_W_m2K = naturalConvectionVerticalPlateH(Tsurf_C, Tinf_C, L_m, pressure_Pa)
Tfilm_K = 0.5 * (Tsurf_C + Tinf_C) + 273.15;
props = airProps(Tfilm_K, pressure_Pa);
beta_1_K = 1 / max(Tfilm_K, 1e-9);
nu_m2_s = props.mu_Pa_s / max(props.rho_kg_m3, 1e-12);
alpha_m2_s = props.k_W_mK / max(props.rho_kg_m3 * props.cp_J_kgK, 1e-12);
g_m_s2 = 9.81;
deltaT_K = max(abs(Tsurf_C - Tinf_C), 1e-6);

% [R4] Ra = g * beta * DeltaT * L^3 / (nu * alpha)
Ra = g_m_s2 * beta_1_K * deltaT_K * L_m^3 / max(nu_m2_s * alpha_m2_s, 1e-12);

% [R4] Churchill-Chu all-Ra vertical plate correlation.
Nu = (0.825 + 0.387 * Ra^(1/6) / ...
    (1 + (0.492 / max(props.Pr, 1e-9))^(9/16))^(8/27))^2;

% [R4] h = Nu*k/L
h_W_m2K = Nu * props.k_W_mK / max(L_m, 1e-9);
end

function h_W_m2K = naturalConvectionHorizontalPlateH(Tsurf_C, Tinf_C, L_m, pressure_Pa)
Tfilm_K = 0.5 * (Tsurf_C + Tinf_C) + 273.15;
props = airProps(Tfilm_K, pressure_Pa);
beta_1_K = 1 / max(Tfilm_K, 1e-9);
nu_m2_s = props.mu_Pa_s / max(props.rho_kg_m3, 1e-12);
alpha_m2_s = props.k_W_mK / max(props.rho_kg_m3 * props.cp_J_kgK, 1e-12);
g_m_s2 = 9.81;
deltaT_K = max(abs(Tsurf_C - Tinf_C), 1e-6);

% [R5] Horizontal plate free convection using standard hot-up/cold-down
% and hot-down/cold-up correlations based on the sign of (Tsurf - Tinf).
Ra = g_m_s2 * beta_1_K * deltaT_K * L_m^3 / max(nu_m2_s * alpha_m2_s, 1e-12);
if Tsurf_C >= Tinf_C
    if Ra <= 1e7
        Nu = 0.54 * Ra^(1/4);
    else
        Nu = 0.15 * Ra^(1/3);
    end
else
    Nu = 0.27 * Ra^(1/4);
end

h_W_m2K = Nu * props.k_W_mK / max(L_m, 1e-9);
end

function flowMin_Lpm = minimumTotalFlowForOutletLimit(params, Qreq_W, dutyFraction)
deltaTuse_K = params.limits.maxAirOutletTemp_C - params.bin.setpoint_C;
if deltaTuse_K <= 0 || dutyFraction <= 0
    flowMin_Lpm = inf;
    return;
end

propsIn = airProps( ...
    params.environment.greenhouseAir_C + 273.15, ...
    params.environment.pressure_Pa);

% [R5] Q = m_dot * cp * (T_injected - T_bed)
mdotMin_kg_s = (Qreq_W / dutyFraction) / max(propsIn.cp_J_kgK * deltaTuse_K, 1e-9);
flowMin_Lpm = mdotMin_kg_s / max(propsIn.rho_kg_m3, 1e-9) * 60 * 1000;
end

function operatingPoints = runOperatingSweep(params)
nPitch = numel(params.sweep.pitch_mm);
nFlow = numel(params.sweep.totalFlow_Lpm);
nVolt = numel(params.sweep.voltage_V);
nPts = nPitch * nFlow * nVolt;
operatingPoints = runOperatingSweepSubset( ...
    params, params.sweep.pitch_mm, params.sweep.totalFlow_Lpm, ...
    params.sweep.voltage_V, nPts);
end

function operatingPoints = runTargetedOperatingSweep(params, targetRec)
if nargin < 2 || isempty(fieldnames(targetRec))
    operatingPoints = runOperatingSweep(params);
    return;
end

pitchVals_mm = nearestSweepValue(params.sweep.pitch_mm, targetRec.pitch_mm);
flowVals_Lpm = wireTradeStudyFlowGrid(params);
voltageCenter_V = wireTradeStudyVoltageCenter(params, targetRec);
voltVals_V = localSweepWindow( ...
    params.sweep.voltage_V, voltageCenter_V, ...
    params.sweep.tradeStudyVoltageHalfWindow);
nPts = numel(pitchVals_mm) * numel(flowVals_Lpm) * numel(voltVals_V);
operatingPoints = runOperatingSweepSubset( ...
    params, pitchVals_mm, flowVals_Lpm, voltVals_V, nPts);
end

function flowVals_Lpm = wireTradeStudyFlowGrid(params)
base = params.sweep.totalFlow_Lpm(:).';
nPts = max(3, round(params.sweep.tradeStudyFlowPointCount));
if numel(base) <= nPts
    flowVals_Lpm = base;
    return;
end
idx = unique(round(linspace(1, numel(base), nPts)));
flowVals_Lpm = base(idx);
end

function voltageCenter_V = wireTradeStudyVoltageCenter(params, targetRec)
voltageCenter_V = targetRec.voltage_V;
if ~isfield(targetRec, 'totalCurrent_A') || ...
        ~isfinite(targetRec.totalCurrent_A) || targetRec.totalCurrent_A <= 0
    return;
end

nTubes = max(1, params.aeration.nParallelTubes);
targetBranchCurrent_A = targetRec.totalCurrent_A / nTubes;
wire = buildWireGeometry(params);
wireArea_m2 = pi * params.wire.diameter_m^2 / 4;
branchResistance20_Ohm = params.wire.rho20_Ohm_m * wire.length_m / ...
    max(wireArea_m2, 1e-12);
trialVoltage_V = targetBranchCurrent_A * branchResistance20_Ohm;
if isfinite(trialVoltage_V) && trialVoltage_V > 0
    voltageCenter_V = trialVoltage_V;
end
end

function operatingPoints = runOperatingSweepSubset( ...
        params, pitchVals_mm, flowVals_Lpm, voltVals_V, nPts)
if nargin < 5
    nPts = numel(pitchVals_mm) * numel(flowVals_Lpm) * numel(voltVals_V);
end

template = struct( ...
    'pitch_mm', NaN, ...
    'nParallelTubes', NaN, ...
    'parallelStringCount', NaN, ...
    'voltage_V', NaN, ...
    'totalFlow_Lpm', NaN, ...
    'branchFlow_Lpm', NaN, ...
    'wireLength_m', NaN, ...
    'branchResistance_Ohm', NaN, ...
    'supplyEquivalentResistance_Ohm', NaN, ...
    'branchCurrent_A', NaN, ...
    'totalCurrent_A', NaN, ...
    'branchPower_W', NaN, ...
    'totalPower_W', NaN, ...
    'hWire_W_m2K', NaN, ...
    'airOutlet_C', NaN, ...
    'wireMax_C', NaN, ...
    'branchQtoBed_W', NaN, ...
    'QtoBed_W', NaN, ...
    'deltaP_Pa', NaN, ...
    'distributionDeltaP_Pa', NaN, ...
    'distributionHeader_Pa', NaN, ...
    'distributionHeaderToSplitter_Pa', NaN, ...
    'distributionSplitterBody_Pa', NaN, ...
    'distributionConnectorFriction_Pa', NaN, ...
    'distributionConnectorToBranch_Pa', NaN, ...
    'flowPerHole_Lpm', NaN, ...
    'holeVelocity_m_s', NaN, ...
    'waterLoss_kg_day', NaN, ...
    'latentEvap_W', NaN, ...
    'electricalTopologyLabel', '', ...
    'stringTubeCounts', [], ...
    'tubeCurrentMin_A', NaN, ...
    'tubeCurrentMean_A', NaN, ...
    'tubeCurrentMax_A', NaN, ...
    'tubePowerMin_W', NaN, ...
    'tubePowerMean_W', NaN, ...
    'tubePowerMax_W', NaN);
operatingPoints = repmat(template, nPts, 1);
nPitch = numel(pitchVals_mm);
nFlow = numel(flowVals_Lpm);
nRows = nPitch * nFlow;
rowResults = cell(nRows, 1);
useParallel = parallelSweepEnabled(params.parallel);

if useParallel
    ensureParallelPool(params.parallel);
    useParallel = ~isempty(gcp('nocreate'));
end

if useParallel
    parfor rowIdx = 1:nRows
        [iPitch, iFlow] = ind2sub([nPitch, nFlow], rowIdx);
        rowResults{rowIdx} = runOperatingSweepRow( ...
            params, pitchVals_mm(iPitch), flowVals_Lpm(iFlow), voltVals_V, template);
    end
else
    for rowIdx = 1:nRows
        [iPitch, iFlow] = ind2sub([nPitch, nFlow], rowIdx);
        rowResults{rowIdx} = runOperatingSweepRow( ...
            params, pitchVals_mm(iPitch), flowVals_Lpm(iFlow), voltVals_V, template);
    end
end

idx = 0;
for rowIdx = 1:nRows
    rowOps = rowResults{rowIdx};
    nRowOps = numel(rowOps);
    operatingPoints(idx + (1:nRowOps)) = rowOps;
    idx = idx + nRowOps;
end
end

function rowOps = runOperatingSweepRow( ...
        params, pitch_mm, totalFlow_Lpm, voltVals_V, template)
trialParams = params;
trialParams.heaterTube.pitch_m = pitch_mm / 1000;
wire = buildWireGeometry(trialParams);
rowOps = repmat(template, numel(voltVals_V), 1);
prevAvgWire_C = NaN;

for iVolt = 1:numel(voltVals_V)
    voltage_V = voltVals_V(iVolt);
    op = evaluateOperatingPoint( ...
        trialParams, voltage_V, totalFlow_Lpm, ...
        params.bin.setpoint_C, params.environment.greenhouseAir_C, ...
        prevAvgWire_C, wire);
    prevAvgWire_C = op.branch.heater.wireMeanTemp_C;

    rowOps(iVolt).pitch_mm = 1000 * wire.pitch_m;
    rowOps(iVolt).nParallelTubes = params.aeration.nParallelTubes;
    rowOps(iVolt).parallelStringCount = op.total.parallelStringCount;
    rowOps(iVolt).voltage_V = voltage_V;
    rowOps(iVolt).totalFlow_Lpm = totalFlow_Lpm;
    rowOps(iVolt).branchFlow_Lpm = op.branch.flow_Lpm;
    rowOps(iVolt).wireLength_m = op.branch.wire.length_m;
    rowOps(iVolt).branchResistance_Ohm = op.branch.wire.resistance_Ohm;
    rowOps(iVolt).supplyEquivalentResistance_Ohm = op.total.equivalentResistance_Ohm;
    rowOps(iVolt).branchCurrent_A = op.branch.wire.current_A;
    rowOps(iVolt).totalCurrent_A = op.total.current_A;
    rowOps(iVolt).branchPower_W = op.branch.wire.power_W;
    rowOps(iVolt).totalPower_W = op.total.power_W;
    rowOps(iVolt).hWire_W_m2K = op.branch.heater.hWireMean_W_m2K;
    rowOps(iVolt).airOutlet_C = op.branch.heater.airOutletTemp_C;
    rowOps(iVolt).wireMax_C = op.branch.heater.wireMaxTemp_C;
    rowOps(iVolt).branchQtoBed_W = op.branch.aeration.QtoBed_W;
    rowOps(iVolt).QtoBed_W = op.total.QtoBed_W;
    rowOps(iVolt).deltaP_Pa = op.total.deltaP_total_Pa;
    rowOps(iVolt).distributionDeltaP_Pa = op.total.distributionDeltaP_Pa;
    rowOps(iVolt).distributionHeader_Pa = op.total.distributionHeader_Pa;
    rowOps(iVolt).distributionHeaderToSplitter_Pa = op.total.distributionHeaderToSplitter_Pa;
    rowOps(iVolt).distributionSplitterBody_Pa = op.total.distributionSplitterBody_Pa;
    rowOps(iVolt).distributionConnectorFriction_Pa = op.total.distributionConnectorFriction_Pa;
    rowOps(iVolt).distributionConnectorToBranch_Pa = op.total.distributionConnectorToBranch_Pa;
    rowOps(iVolt).flowPerHole_Lpm = op.total.flowPerHole_Lpm;
    rowOps(iVolt).holeVelocity_m_s = op.total.holeVelocity_m_s;
    rowOps(iVolt).waterLoss_kg_day = op.total.waterLoss_kg_day;
    rowOps(iVolt).latentEvap_W = op.total.latentEvap_W;
    rowOps(iVolt).electricalTopologyLabel = op.total.electricalTopologyLabel;
    rowOps(iVolt).stringTubeCounts = op.total.stringTubeCounts;
    rowOps(iVolt).tubeCurrentMin_A = op.total.tubeCurrentMin_A;
    rowOps(iVolt).tubeCurrentMean_A = op.total.tubeCurrentMean_A;
    rowOps(iVolt).tubeCurrentMax_A = op.total.tubeCurrentMax_A;
    rowOps(iVolt).tubePowerMin_W = op.total.tubePowerMin_W;
    rowOps(iVolt).tubePowerMean_W = op.total.tubePowerMean_W;
    rowOps(iVolt).tubePowerMax_W = op.total.tubePowerMax_W;
end
end

function tf = parallelSweepEnabled(parallelCfg)
tf = false;
if ~isstruct(parallelCfg) || ~isfield(parallelCfg, 'enabled') || ...
        ~parallelCfg.enabled
    return;
end
try
    tf = license('test', 'Distrib_Computing_Toolbox');
catch
    tf = false;
end
end

function ensureParallelPool(parallelCfg)
if ~isstruct(parallelCfg)
    return;
end
if isfield(parallelCfg, 'autoStartPool') && ~parallelCfg.autoStartPool
    return;
end

pool = gcp('nocreate');
if ~isempty(pool)
    return;
end

workers = 0;
if isfield(parallelCfg, 'workers') && ~isempty(parallelCfg.workers)
    workers = max(0, round(parallelCfg.workers));
end

try
    if workers > 0
        parpool('local', workers);
    else
        parpool('local');
    end
catch
    % Fall back silently; the sweep will still run serially if the pool
    % cannot be started.
end
end

function value = nearestSweepValue(values, target)
[~, idx] = min(abs(values - target));
value = values(idx);
end

function values = localSweepWindow(gridValues, targetValue, halfWindow)
[~, centerIdx] = min(abs(gridValues - targetValue));
lo = max(1, centerIdx - max(0, halfWindow));
hi = min(numel(gridValues), centerIdx + max(0, halfWindow));
values = gridValues(lo:hi);
end

function designPoints = applyScenarioToSweep(operatingPoints, binModel, Qreq_W, params)
designPoints = operatingPoints;

for i = 1:numel(designPoints)
    [TbottomFull_C, TtopFull_C] = twoNodeSteadyState( ...
        binModel, params.environment.greenhouseAir_C, ...
        designPoints(i).QtoBed_W, params.bin.heatToBottomFraction);

    spread_C = abs(TbottomFull_C - TtopFull_C);
    dutyNeeded = Qreq_W / max(designPoints(i).QtoBed_W, 1e-9);
    constraintSafe = ...
        designPoints(i).totalCurrent_A <= params.limits.maxTotalCurrent_A && ...
        designPoints(i).wireMax_C <= params.limits.maxWireTemp_C && ...
        designPoints(i).airOutlet_C <= params.limits.maxAirOutletTemp_C && ...
        designPoints(i).deltaP_Pa <= params.limits.maxPressureDrop_Pa && ...
        designPoints(i).waterLoss_kg_day <= params.limits.maxWaterLoss_kg_day && ...
        designPoints(i).holeVelocity_m_s <= params.limits.maxHoleVelocity_m_s && ...
        TbottomFull_C >= params.optimization.minBottomTemp_C && ...
        TtopFull_C >= params.optimization.minTopTemp_C && ...
        spread_C <= params.optimization.spreadLimit_C && ...
        designPoints(i).QtoBed_W > 0;
    feasible = constraintSafe && dutyNeeded <= 1.0;

    designPoints(i).dutyNeeded = dutyNeeded;
    designPoints(i).TbottomFullPower_C = TbottomFull_C;
    designPoints(i).TtopFullPower_C = TtopFull_C;
    designPoints(i).isConstraintSafe = constraintSafe;
    designPoints(i).isFeasible = feasible;
end
end

function recommended = selectRecommendedPoint(designPoints, params)
recommended = struct();
if isempty(designPoints)
    return;
end

constraintSafe = [designPoints.isConstraintSafe];
if ~any(constraintSafe)
    return;
end

[recommended, ~] = selectBestPointByCriteria(designPoints(constraintSafe), params);
recommended.meetsCurrentCap = ...
    recommended.totalCurrent_A <= params.limits.maxTotalCurrent_A;
end

function bestAvailable = selectBestAvailablePoint(designPoints, params)
bestAvailable = struct();
if isempty(designPoints)
    return;
end

[bestAvailable, ~] = selectBestPointByCriteria(designPoints, params);
bestAvailable.meetsCurrentCap = ...
    bestAvailable.totalCurrent_A <= params.limits.maxTotalCurrent_A;
end

function op = evaluateOperatingPoint( ...
        params, voltage_V, totalFlow_Lpm, Tbed_C, Tenv_C, initAvgWire_C, wire)
nTubes = max(1, params.aeration.nParallelTubes);
stringCounts = params.electrical.stringTubeCounts;
if isempty(stringCounts) || sum(stringCounts) ~= nTubes
    stringCounts = preferredElectricalStringTubeCounts(nTubes, params.electrical);
end
branchFlow_Lpm = totalFlow_Lpm / nTubes;

airIn_C = Tenv_C;
airIn_K = airIn_C + 273.15;

if nargin < 6 || isempty(initAvgWire_C) || ~isfinite(initAvgWire_C)
    avgWire_C = airIn_C + 20;
else
    avgWire_C = initAvgWire_C;
end

if nargin < 7 || isempty(wire)
    wire = buildWireGeometry(params);
end

if ~isfield(wire, 'isValid') || ~wire.isValid
    op = invalidOperatingPoint(params, voltage_V, totalFlow_Lpm, wire);
    return;
end

for iter = 1:8
    wireProps = nichrome80Props(params.wire, avgWire_C);
    tubeResistance_Ohm = wireProps.resistance_Ohm_m * wire.length_m;
    seriesCounts = unique(stringCounts);
    groupStates = repmat(struct(), numel(seriesCounts), 1);
    totalWireMean_C = 0;
    totalTubeWeight = 0;

    for iGroup = 1:numel(seriesCounts)
        seriesCount = seriesCounts(iGroup);
        nStringsGroup = sum(stringCounts == seriesCount);
        nTubesGroup = nStringsGroup * seriesCount;
        tubeCurrent_A = voltage_V / max(seriesCount * tubeResistance_Ohm, 1e-9);
        tubePower_W = tubeCurrent_A^2 * tubeResistance_Ohm;
        heater = solveHeaterTube( ...
            params, wire, tubePower_W, branchFlow_Lpm, airIn_K, Tenv_C);

        groupStates(iGroup).seriesCount = seriesCount;
        groupStates(iGroup).nStrings = nStringsGroup;
        groupStates(iGroup).nTubes = nTubesGroup;
        groupStates(iGroup).tubeCurrent_A = tubeCurrent_A;
        groupStates(iGroup).tubePower_W = tubePower_W;
        groupStates(iGroup).heater = heater;

        totalWireMean_C = totalWireMean_C + nTubesGroup * mean(heater.wireProfile_C);
        totalTubeWeight = totalTubeWeight + nTubesGroup;
    end

    newAvgWire_C = totalWireMean_C / max(totalTubeWeight, 1);
    if abs(newAvgWire_C - avgWire_C) < 0.05
        break;
    end
    avgWire_C = 0.65 * avgWire_C + 0.35 * newAvgWire_C;
end

propsAtInlet = airProps(airIn_K, params.environment.pressure_Pa);
volFlowBranch_m3_s = branchFlow_Lpm / (1000 * 60);
mdotBranch_kg_s = volFlowBranch_m3_s * propsAtInlet.rho_kg_m3;
distribution = distributionNetworkLosses( ...
    params, totalFlow_Lpm, branchFlow_Lpm, airIn_K, params.heaterTube.ID_m);

heaterWireProfileSum_C = [];
heaterAirProfileSum_C = [];
heaterHWireWeighted = 0;
heaterWireMeanWeighted = 0;
heaterDeltaPMax_Pa = -Inf;
airOutletMax_C = -Inf;
wireMax_C = -Inf;
branchDeltaP_Pa = -Inf;
flowPerHole_Lpm = -Inf;
holeVelocity_m_s = -Inf;
totalPower_W = 0;
totalCurrent_A = 0;
QtoBedTotal_W = 0;
waterLossTotal_kg_day = 0;
latentEvapTotal_W = 0;
tubeCurrentMin_A = Inf;
tubeCurrentMax_A = -Inf;
tubeCurrentWeighted_A = 0;
tubePowerMin_W = Inf;
tubePowerMax_W = -Inf;
tubePowerWeighted_W = 0;
aerationQtoBedWeighted_W = 0;
aerationDeltaPMax_Pa = -Inf;

for iGroup = 1:numel(groupStates)
    heater = groupStates(iGroup).heater;
    seriesCount = groupStates(iGroup).seriesCount;
    nStringsGroup = groupStates(iGroup).nStrings;
    nTubesGroup = groupStates(iGroup).nTubes;
    tubeCurrent_A = groupStates(iGroup).tubeCurrent_A;
    tubePower_W = groupStates(iGroup).tubePower_W;

    branchAeration = solveAerationTube( ...
        params, heater.airOutletTemp_C, Tbed_C, mdotBranch_kg_s);
    dpExpansion_Pa = abruptExpansionLoss( ...
        mdotBranch_kg_s, heater.airOutletTemp_C + 273.15, ...
        params.heaterTube.ID_m, params.aeration.ID_m, params.environment.pressure_Pa);
    groupBranchDeltaP_Pa = heater.deltaP_Pa + dpExpansion_Pa + branchAeration.deltaP_Pa;

    if isempty(heaterWireProfileSum_C)
        heaterWireProfileSum_C = zeros(size(heater.wireProfile_C));
        heaterAirProfileSum_C = zeros(size(heater.airProfile_C));
    end
    heaterWireProfileSum_C = heaterWireProfileSum_C + nTubesGroup * heater.wireProfile_C;
    heaterAirProfileSum_C = heaterAirProfileSum_C + nTubesGroup * heater.airProfile_C;

    heaterHWireWeighted = heaterHWireWeighted + nTubesGroup * heater.hWireMean_W_m2K;
    heaterWireMeanWeighted = heaterWireMeanWeighted + nTubesGroup * heater.wireMeanTemp_C;
    heaterDeltaPMax_Pa = max(heaterDeltaPMax_Pa, heater.deltaP_Pa);
    airOutletMax_C = max(airOutletMax_C, heater.airOutletTemp_C);
    wireMax_C = max(wireMax_C, heater.wireMaxTemp_C);
    branchDeltaP_Pa = max(branchDeltaP_Pa, groupBranchDeltaP_Pa);
    aerationDeltaPMax_Pa = max(aerationDeltaPMax_Pa, branchAeration.deltaP_Pa);
    flowPerHole_Lpm = max(flowPerHole_Lpm, branchAeration.flowPerHole_Lpm);
    holeVelocity_m_s = max(holeVelocity_m_s, branchAeration.holeVelocity_m_s);

    totalCurrent_A = totalCurrent_A + nStringsGroup * tubeCurrent_A;
    totalPower_W = totalPower_W + nTubesGroup * tubePower_W;
    QtoBedTotal_W = QtoBedTotal_W + nTubesGroup * branchAeration.QtoBed_W;
    waterLossTotal_kg_day = waterLossTotal_kg_day + nTubesGroup * branchAeration.waterLoss_kg_day;
    latentEvapTotal_W = latentEvapTotal_W + nTubesGroup * branchAeration.latentEvap_W;
    aerationQtoBedWeighted_W = aerationQtoBedWeighted_W + nTubesGroup * branchAeration.QtoBed_W;

    tubeCurrentMin_A = min(tubeCurrentMin_A, tubeCurrent_A);
    tubeCurrentMax_A = max(tubeCurrentMax_A, tubeCurrent_A);
    tubeCurrentWeighted_A = tubeCurrentWeighted_A + nTubesGroup * tubeCurrent_A;
    tubePowerMin_W = min(tubePowerMin_W, tubePower_W);
    tubePowerMax_W = max(tubePowerMax_W, tubePower_W);
    tubePowerWeighted_W = tubePowerWeighted_W + nTubesGroup * tubePower_W;
end

meanTubeCurrent_A = tubeCurrentWeighted_A / max(nTubes, 1);
meanTubePower_W = tubePowerWeighted_W / max(nTubes, 1);
meanHeaterHWire_W_m2K = heaterHWireWeighted / max(nTubes, 1);
meanWireTemp_C = heaterWireMeanWeighted / max(nTubes, 1);
meanBranchQtoBed_W = aerationQtoBedWeighted_W / max(nTubes, 1);
equivalentResistance_Ohm = voltage_V / max(totalCurrent_A, 1e-9);
topologyLabel = electricalTopologyLabel(stringCounts);

op = struct();
op.branch = struct();
op.branch.flow_Lpm = branchFlow_Lpm;
op.branch.wire = struct();
op.branch.wire.length_m = wire.length_m;
op.branch.wire.area_m2 = wire.surfaceArea_m2;
op.branch.wire.resistance_Ohm = tubeResistance_Ohm;
op.branch.wire.current_A = meanTubeCurrent_A;
op.branch.wire.currentMin_A = tubeCurrentMin_A;
op.branch.wire.currentMax_A = tubeCurrentMax_A;
op.branch.wire.power_W = meanTubePower_W;
op.branch.wire.powerMin_W = tubePowerMin_W;
op.branch.wire.powerMax_W = tubePowerMax_W;
op.branch.heater = struct();
op.branch.heater.hWireMean_W_m2K = meanHeaterHWire_W_m2K;
op.branch.heater.airOutletTemp_C = airOutletMax_C;
op.branch.heater.wireMaxTemp_C = wireMax_C;
op.branch.heater.wireMeanTemp_C = meanWireTemp_C;
op.branch.heater.wireProfile_C = heaterWireProfileSum_C / max(nTubes, 1);
op.branch.heater.airProfile_C = heaterAirProfileSum_C / max(nTubes, 1);
op.branch.heater.deltaP_Pa = heaterDeltaPMax_Pa;
op.branch.aeration = struct();
op.branch.aeration.QtoBed_W = meanBranchQtoBed_W;
op.branch.aeration.deltaP_Pa = aerationDeltaPMax_Pa;
op.branch.aeration.flowPerHole_Lpm = flowPerHole_Lpm;
op.branch.aeration.holeVelocity_m_s = holeVelocity_m_s;
op.branch.aeration.waterLoss_kg_day = waterLossTotal_kg_day / max(nTubes, 1);
op.branch.aeration.latentEvap_W = latentEvapTotal_W / max(nTubes, 1);
op.branch.deltaP_Pa = branchDeltaP_Pa;

op.total = struct();
op.total.nParallelTubes = nTubes;
op.total.parallelStringCount = numel(stringCounts);
op.total.stringTubeCounts = stringCounts;
op.total.electricalTopologyLabel = topologyLabel;
op.total.flow_Lpm = totalFlow_Lpm;
op.total.mdot_kg_s = nTubes * mdotBranch_kg_s;
op.total.equivalentResistance_Ohm = equivalentResistance_Ohm;
op.total.current_A = totalCurrent_A;
op.total.power_W = totalPower_W;
op.total.QtoBed_W = QtoBedTotal_W;
op.total.distributionDeltaP_Pa = distribution.total_Pa;
op.total.distributionHeader_Pa = distribution.header_Pa;
op.total.distributionHeaderToSplitter_Pa = distribution.headerToSplitter_Pa;
op.total.distributionSplitterBody_Pa = distribution.splitterBody_Pa;
op.total.distributionConnectorFriction_Pa = distribution.connectorFriction_Pa;
op.total.distributionConnectorToBranch_Pa = distribution.connectorToBranch_Pa;
op.total.deltaP_total_Pa = branchDeltaP_Pa + distribution.total_Pa;
op.total.flowPerHole_Lpm = flowPerHole_Lpm;
op.total.holeVelocity_m_s = holeVelocity_m_s;
op.total.waterLoss_kg_day = waterLossTotal_kg_day;
op.total.latentEvap_W = latentEvapTotal_W;
op.total.tubeCurrentMin_A = tubeCurrentMin_A;
op.total.tubeCurrentMean_A = meanTubeCurrent_A;
op.total.tubeCurrentMax_A = tubeCurrentMax_A;
op.total.tubePowerMin_W = tubePowerMin_W;
op.total.tubePowerMean_W = meanTubePower_W;
op.total.tubePowerMax_W = tubePowerMax_W;
end

function wire = buildWireGeometry(params)
coilSpan_m = params.heaterTube.coilSpan_m;
minPitch_m = params.wire.diameter_m * 1.05;
manualWireLength_m = params.heaterTube.manualWireLength_m;
coilMeanD_m = params.heaterTube.coilMeanD_m;
fitClearance_m = params.heaterTube.ID_m - ...
    (coilMeanD_m + params.wire.diameter_m);

wire = struct();
wire.isValid = coilSpan_m > 0 && coilMeanD_m > 0 && fitClearance_m >= 0;
wire.fitClearance_m = fitClearance_m;
wire.outerDiameter_m = coilMeanD_m + params.wire.diameter_m;
wire.turns = NaN;
wire.pitch_m = NaN;
wire.length_m = NaN;
wire.crossSection_m2 = pi * params.wire.diameter_m^2 / 4;
wire.surfaceArea_m2 = NaN;
wire.surfaceAreaPerAxialLength_m = NaN;
wire.turnLength_m = NaN;
if ~wire.isValid
    return;
end

if isfinite(manualWireLength_m) && manualWireLength_m > coilSpan_m
    pitch_m = impliedPitchForWireLength( ...
        manualWireLength_m, coilSpan_m, coilMeanD_m);
    pitch_m = max(pitch_m, minPitch_m);
else
    pitch_m = max(params.heaterTube.pitch_m, minPitch_m);
end
turns = max(1, coilSpan_m / pitch_m);

% [R7] Helix geometry: L_turn = sqrt((pi*D_mean)^2 + pitch^2)
turnLength_m = hypot(pi * coilMeanD_m, pitch_m);
wireLength_m = turns * turnLength_m;

wire.turns = turns;
wire.pitch_m = pitch_m;
wire.length_m = wireLength_m;
wire.surfaceArea_m2 = pi * params.wire.diameter_m * wireLength_m;
wire.surfaceAreaPerAxialLength_m = wire.surfaceArea_m2 / max(coilSpan_m, 1e-9);
wire.turnLength_m = turnLength_m;
end

function op = invalidOperatingPoint(params, voltage_V, totalFlow_Lpm, wire)
nTubes = max(1, params.aeration.nParallelTubes);
stringCounts = params.electrical.stringTubeCounts;
if isempty(stringCounts) || sum(stringCounts) ~= nTubes
    stringCounts = preferredElectricalStringTubeCounts(nTubes, params.electrical);
end
branchFlow_Lpm = totalFlow_Lpm / nTubes;

heater = struct();
heater.hWireMean_W_m2K = NaN;
heater.airOutletTemp_C = params.environment.greenhouseAir_C;
heater.wireMaxTemp_C = Inf;
heater.wireMeanTemp_C = NaN;
heater.wireProfile_C = nan(max(params.heaterTube.segments, 1), 1);
heater.airProfile_C = nan(max(params.heaterTube.segments + 1, 1), 1);
heater.deltaP_Pa = Inf;

aeration = struct();
aeration.QtoBed_W = 0;
aeration.deltaP_Pa = Inf;
aeration.flowPerHole_Lpm = NaN;
aeration.holeVelocity_m_s = Inf;
aeration.waterLoss_kg_day = Inf;
aeration.latentEvap_W = NaN;

op = struct();
op.branch = struct();
op.branch.flow_Lpm = branchFlow_Lpm;
op.branch.wire = struct();
op.branch.wire.length_m = wire.length_m;
op.branch.wire.area_m2 = wire.surfaceArea_m2;
op.branch.wire.resistance_Ohm = Inf;
op.branch.wire.current_A = Inf;
op.branch.wire.power_W = Inf;
op.branch.heater = heater;
op.branch.aeration = aeration;
op.branch.deltaP_Pa = Inf;

op.total = struct();
op.total.nParallelTubes = nTubes;
op.total.parallelStringCount = numel(stringCounts);
op.total.stringTubeCounts = stringCounts;
op.total.electricalTopologyLabel = electricalTopologyLabel(stringCounts);
op.total.flow_Lpm = totalFlow_Lpm;
op.total.mdot_kg_s = NaN;
op.total.equivalentResistance_Ohm = Inf;
op.total.current_A = Inf;
op.total.power_W = Inf;
op.total.QtoBed_W = 0;
op.total.distributionDeltaP_Pa = Inf;
op.total.distributionHeader_Pa = Inf;
op.total.distributionHeaderToSplitter_Pa = Inf;
op.total.distributionSplitterBody_Pa = Inf;
op.total.distributionConnectorFriction_Pa = Inf;
op.total.distributionConnectorToBranch_Pa = Inf;
op.total.deltaP_total_Pa = Inf;
op.total.flowPerHole_Lpm = NaN;
op.total.holeVelocity_m_s = Inf;
op.total.waterLoss_kg_day = Inf;
op.total.latentEvap_W = NaN;
op.total.tubeCurrentMin_A = Inf;
op.total.tubeCurrentMean_A = Inf;
op.total.tubeCurrentMax_A = Inf;
op.total.tubePowerMin_W = Inf;
op.total.tubePowerMean_W = Inf;
op.total.tubePowerMax_W = Inf;
end

function pitch_m = impliedPitchForWireLength( ...
        wireLength_m, coilSpan_m, coilMeanD_m)
circumference_m = pi * coilMeanD_m;
denom = max(wireLength_m^2 - coilSpan_m^2, 1e-12);
pitch_m = coilSpan_m * circumference_m / sqrt(denom);
end

function heater = solveHeaterTube(params, wire, power_W, flow_Lpm, airIn_K, Tamb_C)
N = params.heaterTube.segments;
L_m = params.heaterTube.length_m;
dx_m = L_m / N;
ID_m = params.heaterTube.ID_m;
OD_m = params.heaterTube.OD_m;
Atube_m2 = pi * ID_m^2 / 4;
Poutside_m = pi * OD_m;
Tair_K = airIn_K;

propsIn = airProps(airIn_K, params.environment.pressure_Pa);
volFlow_m3_s = flow_Lpm / (1000 * 60);
mdot_kg_s = volFlow_m3_s * propsIn.rho_kg_m3;

wireProfile_C = zeros(N, 1);
airProfile_C = zeros(N + 1, 1);
airProfile_C(1) = airIn_K - 273.15;
hWireProfile = zeros(N, 1);
dpFriction_Pa = 0;

for i = 1:N
    props = airProps(Tair_K, params.environment.pressure_Pa);

    % [R5] u = m_dot / (rho * A)
    uBulk_m_s = mdot_kg_s / max(props.rho_kg_m3 * Atube_m2, 1e-12);

    helixAngleFactor = (pi * params.heaterTube.coilMeanD_m) / ...
        hypot(pi * params.heaterTube.coilMeanD_m, params.heaterTube.pitch_m);
    uNormal_m_s = max(1e-6, uBulk_m_s * helixAngleFactor);

    % [R1] Re and Nu for crossflow over the wire approximated as a cylinder.
    ReWire = props.rho_kg_m3 * uNormal_m_s * params.wire.diameter_m / ...
        max(props.mu_Pa_s, 1e-12);
    NuWire = churchillBernsteinNu(ReWire, props.Pr);

    % [R1] h = Nu*k/D
    hWire_W_m2K = params.calibration.hWireMultiplier * ...
        NuWire * props.k_W_mK / params.wire.diameter_m;

    % [R2, R5] Internal-tube Re and Nu for the tube wall loss estimate.
    ReTube = props.rho_kg_m3 * uBulk_m_s * ID_m / max(props.mu_Pa_s, 1e-12);
    NuTube = internalTubeNu(ReTube, props.Pr, ID_m, L_m);
    hInsideWall_W_m2K = params.calibration.hWallMultiplier * ...
        NuTube * props.k_W_mK / ID_m;

    Uloss_W_m2K = heaterTubeOverallU(params, hInsideWall_W_m2K);

    % [R5] q_loss = U*A*DeltaT from the heater shell to surroundings.
    qLoss_W = Uloss_W_m2K * Poutside_m * dx_m * ...
        max((Tair_K - 273.15) - Tamb_C, 0);

    qGen_W = power_W * dx_m / max(L_m, 1e-9);
    qToAir_W = qGen_W - qLoss_W;

    % [R5] m_dot*cp*dT = q_to_air
    Tair_K = Tair_K + qToAir_W / max(mdot_kg_s * props.cp_J_kgK, 1e-9);

    % [R5] Local wire-to-fluid temperature rise from Newton's law of cooling.
    deltaTwire_K = qGen_W / max(hWire_W_m2K * wire.surfaceAreaPerAxialLength_m * dx_m, 1e-9);
    wireProfile_C(i) = (Tair_K - 273.15) + deltaTwire_K;
    airProfile_C(i + 1) = Tair_K - 273.15;
    hWireProfile(i) = hWire_W_m2K;

    % [R3] Darcy-Weisbach pressure drop with Churchill explicit f.
    f = churchillFrictionFactor(ReTube);
    dpFriction_Pa = dpFriction_Pa + ...
        (f * dx_m / ID_m + params.heaterTube.extraMinorK / N) * ...
        0.5 * props.rho_kg_m3 * uBulk_m_s^2;
end

heater = struct();
heater.airInletTemp_C = airIn_K - 273.15;
heater.airOutletTemp_C = Tair_K - 273.15;
heater.airProfile_C = airProfile_C;
heater.wireProfile_C = wireProfile_C;
heater.wireMaxTemp_C = max(wireProfile_C);
heater.wireMeanTemp_C = mean(wireProfile_C);
heater.hWireMean_W_m2K = mean(hWireProfile);
heater.deltaP_Pa = dpFriction_Pa;
end

function U_W_m2K = heaterTubeOverallU(params, hInsideWall_W_m2K)
% [R5] Radial resistance network approximated as series resistances.
Rin = 1 / max(hInsideWall_W_m2K, 1e-9);
Rwall = (params.heaterTube.OD_m - params.heaterTube.ID_m) / ...
    (2 * max(params.heaterTube.wallK_W_mK, 1e-9));
Rins = params.heaterTube.insulationThickness_m / ...
    max(params.heaterTube.insulationK_W_mK, 1e-9);
Rout = 1 / max(params.heaterTube.externalH_W_m2K, 1e-9);
U_W_m2K = 1 / (Rin + Rwall + Rins + Rout);
end

function aeration = solveAerationTube(params, TairIn_C, Tbed_C, mdotIn_kg_s)
propsRef = airProps(TairIn_C + 273.15, params.environment.pressure_Pa);
totalHoleArea_m2 = holesPerTube(params) * perforationHoleArea(params);
Cd = perforationDischargeCoefficient(params);
massBalanceTol_kg_s = max(1e-9, 1e-4 * mdotIn_kg_s);

if mdotIn_kg_s <= 0 || totalHoleArea_m2 <= 0 || Cd <= 0
    aeration = marchAerationTubeClosedEnd(params, TairIn_C, Tbed_C, mdotIn_kg_s, 0.0);
    aeration.deltaP_Pa = 0.0;
    aeration.inletGaugePressure_Pa = 0.0;
    return
end

hi_Pa = max((mdotIn_kg_s / max(Cd * totalHoleArea_m2, 1e-12))^2 / ...
    max(2 * max(propsRef.rho_kg_m3, 1e-12), 1e-12), 1.0);
stateHi = marchAerationTubeClosedEnd(params, TairIn_C, Tbed_C, mdotIn_kg_s, hi_Pa);

for iter = 1:40
    if stateHi.remainingMdot_kg_s <= massBalanceTol_kg_s
        break
    end
    hi_Pa = 2 * hi_Pa;
    stateHi = marchAerationTubeClosedEnd(params, TairIn_C, Tbed_C, mdotIn_kg_s, hi_Pa);
end

if stateHi.remainingMdot_kg_s > massBalanceTol_kg_s
    aeration = stateHi;
    aeration.deltaP_Pa = hi_Pa;
    aeration.inletGaugePressure_Pa = hi_Pa;
    return
end

lo_Pa = 0.0;
bestPressure_Pa = hi_Pa;
bestState = stateHi;
for iter = 1:80
    mid_Pa = 0.5 * (lo_Pa + hi_Pa);
    stateMid = marchAerationTubeClosedEnd(params, TairIn_C, Tbed_C, mdotIn_kg_s, mid_Pa);
    if stateMid.remainingMdot_kg_s <= massBalanceTol_kg_s
        hi_Pa = mid_Pa;
        bestPressure_Pa = mid_Pa;
        bestState = stateMid;
    else
        lo_Pa = mid_Pa;
    end
    if hi_Pa - lo_Pa <= max(1e-6, 1e-6 * max(bestPressure_Pa, 1.0))
        break
    end
end

aeration = bestState;
aeration.deltaP_Pa = bestPressure_Pa;
aeration.inletGaugePressure_Pa = bestPressure_Pa;
end

function aeration = marchAerationTubeClosedEnd(params, TairIn_C, Tbed_C, mdotIn_kg_s, inletGaugePressure_Pa)
N = params.aeration.segments;
L_m = params.aeration.length_m;
dx_m = L_m / max(N, 1);
ID_m = params.aeration.ID_m;
OD_m = params.aeration.OD_m;
Atube_m2 = pi * ID_m^2 / 4;
Pout_m = pi * OD_m;
holesPerSeg = max(holesPerTube(params) / max(N, 1), 1e-12);
segHoleArea_m2 = holesPerSeg * perforationHoleArea(params);
bedPressureModel = buildSegmentwiseBedPressureModel(params, L_m, N, OD_m);
boundaryProps = airProps(Tbed_C + 273.15, params.environment.pressure_Pa);
bedBoundaryGauge_Pa = derivedOutletBoundaryGaugePressure_Pa( ...
    mdotIn_kg_s, boundaryProps.rho_kg_m3, bedPressureModel);

remainingMdot_kg_s = mdotIn_kg_s;
Tair_C = TairIn_C;
QtoBed_W = 0;
latentEvap_W = 0;
waterLoss_kg_s = 0;
dpTube_Pa = 0;
peakFlowPerHole_Lpm = 0;
peakHoleVelocity_m_s = 0;
localGaugePressure_Pa = inletGaugePressure_Pa;

airProfile_C = zeros(N + 1, 1);
airProfile_C(1) = Tair_C;
pressureProfileGauge_Pa = zeros(N + 1, 1);
pressureProfileGauge_Pa(1) = inletGaugePressure_Pa;
bedPressureProfileGauge_Pa = zeros(N, 1);

for i = 1:N
    localPressure_Pa = params.environment.pressure_Pa + ...
        max(localGaugePressure_Pa, bedBoundaryGauge_Pa);
    props = airProps(Tair_C + 273.15, localPressure_Pa);
    uBulk_m_s = remainingMdot_kg_s / max(props.rho_kg_m3 * Atube_m2, 1e-12);
    ReTube = props.rho_kg_m3 * uBulk_m_s * ID_m / max(props.mu_Pa_s, 1e-12);
    NuTube = internalTubeNu(ReTube, props.Pr, ID_m, L_m);
    hInside_W_m2K = params.calibration.hWallMultiplier * ...
        NuTube * props.k_W_mK / ID_m;
    hOutside_W_m2K = params.calibration.bedHTMultiplier * params.aeration.bedH_W_m2K;

    U_W_m2K = 1 / ( ...
        1 / max(hInside_W_m2K, 1e-9) + ...
        (OD_m - ID_m) / (2 * max(params.aeration.wallK_W_mK, 1e-9)) + ...
        1 / max(hOutside_W_m2K, 1e-9));

    [dm_kg_s, localBedGauge_Pa] = solveSegmentReleaseWithBedBackpressure( ...
        params, localGaugePressure_Pa, props.rho_kg_m3, props.mu_Pa_s, ...
        segHoleArea_m2, remainingMdot_kg_s, bedPressureModel, i, bedBoundaryGauge_Pa);
    bedPressureProfileGauge_Pa(i) = localBedGauge_Pa;
    remainingAfterRelease_kg_s = max(remainingMdot_kg_s - dm_kg_s, 1e-9);
    mdotMean_kg_s = max(remainingMdot_kg_s - 0.5 * dm_kg_s, 1e-12);
    UAseg_W_K = U_W_m2K * Pout_m * dx_m;
    wallFactor = exp(-UAseg_W_K / max(mdotMean_kg_s * props.cp_J_kgK, 1e-12));
    TairAfterWall_C = Tbed_C + (Tair_C - Tbed_C) * wallFactor;
    qWall_W = mdotMean_kg_s * props.cp_J_kgK * (Tair_C - TairAfterWall_C);
    Trelease_C = 0.5 * (Tair_C + TairAfterWall_C);
    qJet_W = dm_kg_s * props.cp_J_kgK * (Trelease_C - Tbed_C);

    releaseProps = airProps(Trelease_C + 273.15, localPressure_Pa);
    volRelease_m3_s = dm_kg_s / max(releaseProps.rho_kg_m3, 1e-12);
    localFlowPerHole_Lpm = volRelease_m3_s * 60 * 1000 / max(holesPerSeg, 1e-12);
    localHoleVelocity_m_s = volRelease_m3_s / max(segHoleArea_m2, 1e-12);
    peakFlowPerHole_Lpm = max(peakFlowPerHole_Lpm, localFlowPerHole_Lpm);
    peakHoleVelocity_m_s = max(peakHoleVelocity_m_s, localHoleVelocity_m_s);

    qAvail_W = max(qJet_W, 0);
    [mEvap_kg_s, qEvap_W] = evaporationLossSegment( ...
        params, dm_kg_s, Trelease_C, Tbed_C, qAvail_W, localHoleVelocity_m_s);
    QtoBed_W = QtoBed_W + qWall_W + qJet_W - qEvap_W;
    latentEvap_W = latentEvap_W + qEvap_W;
    waterLoss_kg_s = waterLoss_kg_s + mEvap_kg_s;

    Tair_C = TairAfterWall_C;
    remainingMdot_kg_s = remainingAfterRelease_kg_s;
    airProfile_C(i + 1) = Tair_C;

    f = churchillFrictionFactor(ReTube);
    segDp_Pa = f * dx_m / ID_m * 0.5 * props.rho_kg_m3 * uBulk_m_s^2;
    dpTube_Pa = dpTube_Pa + segDp_Pa;
    localGaugePressure_Pa = max(localGaugePressure_Pa - segDp_Pa, localBedGauge_Pa);
    pressureProfileGauge_Pa(i + 1) = localGaugePressure_Pa;
end

aeration = struct();
aeration.QtoBed_W = QtoBed_W;
aeration.airOutlet_C = Tair_C;
aeration.airProfile_C = airProfile_C;
aeration.deltaP_Pa = inletGaugePressure_Pa;
aeration.flowPerHole_Lpm = peakFlowPerHole_Lpm;
aeration.holeVelocity_m_s = peakHoleVelocity_m_s;
aeration.waterLoss_kg_day = waterLoss_kg_s * 86400;
aeration.latentEvap_W = latentEvap_W;
aeration.remainingMdot_kg_s = remainingMdot_kg_s;
aeration.tubeFrictionDrop_Pa = dpTube_Pa;
aeration.bedBoundaryGaugePressure_Pa = bedBoundaryGauge_Pa;
aeration.pressureProfileGauge_Pa = pressureProfileGauge_Pa;
aeration.bedPressureProfileGauge_Pa = bedPressureProfileGauge_Pa;
aeration.maxBedGaugePressure_Pa = max(bedPressureProfileGauge_Pa);
end

function holes = totalPerforationHoles(params)
holes = max(1, round(params.aeration.nParallelTubes * params.perforation.holesPerTube));
end

function holes = holesPerTube(params)
holes = max(1, round(params.perforation.holesPerTube));
end

function area_m2 = perforationHoleArea(params)
area_m2 = pi * max(params.perforation.holeDiameter_m, 1e-9)^2 / 4;
end

function Cd = perforationDischargeCoefficient(params)
Cd = 0.60;
if isfield(params, 'perforation') && isstruct(params.perforation) && ...
        isfield(params.perforation, 'dischargeCoefficient') && ...
        ~isempty(params.perforation.dischargeCoefficient)
    Cd = max(params.perforation.dischargeCoefficient, 0.0);
end
end

function bedSideGauge_Pa = aerationBedSideGaugePressure(params)
bedSideGauge_Pa = 0.0;
if isfield(params, 'aeration') && isstruct(params.aeration) && ...
        isfield(params.aeration, 'bedSideGaugePressure_Pa') && ...
        ~isempty(params.aeration.bedSideGaugePressure_Pa)
    bedSideGauge_Pa = params.aeration.bedSideGaugePressure_Pa;
end
end

function centerHeight_m = defaultRepresentativeTubeCenterHeight_m(fillHeight_m, tubeOD_m)
tubeRadius_m = max(0.5 * max(tubeOD_m, 0.0), 0.0);
centerHeight_m = min(max(0.10 * max(fillHeight_m, 0.0), tubeRadius_m), ...
    max(fillHeight_m - tubeRadius_m, tubeRadius_m));
end

function bedModel = buildSegmentwiseBedPressureModel(params, length_m, N, tubeOD_m)
bedModel = struct();
bedModel.enabled = false;
bedModel.mode = 'configured_boundary_only';
bedModel.boundaryGaugePressure_Pa = aerationBedSideGaugePressure(params);
bedModel.boundaryGaugeOffset_Pa = aerationBedSideGaugePressure(params);
bedModel.derivedBoundaryEnabled = false;
bedModel.outletAreaPerBranch_m2 = NaN;
bedModel.outletDischargeCoefficient = NaN;
bedModel.porosity = 0.80;
bedModel.particleDiameter_m = 0.015;
bedModel.particleDiameterSource = 'fallback_constant';
bedModel.representativeTubeCenterHeight_m = NaN;
bedModel.verticalEscapeDistance_m = NaN;
bedModel.tributaryArea_m2 = NaN;
bedModel.segmentPathLengths_m = zeros(max(N, 1), 1);
bedModel.segmentCenters_m = zeros(max(N, 1), 1);

if ~(isfield(params, 'bin') && isstruct(params.bin) && isfield(params, 'aeration') && isstruct(params.aeration))
    return;
end

if isfield(params.aeration, 'bedPorosity') && ~isempty(params.aeration.bedPorosity)
    bedModel.porosity = params.aeration.bedPorosity;
end
particleDiameter_m = NaN;
particleDiameterIsOverride = false;
if isfield(params.aeration, 'bedParticleDiameter_m') && ~isempty(params.aeration.bedParticleDiameter_m)
    particleDiameter_m = params.aeration.bedParticleDiameter_m;
    particleDiameterIsOverride = isfinite(particleDiameter_m) && particleDiameter_m > 0.0;
end
if particleDiameterIsOverride
    bedModel.particleDiameterSource = 'aeration_override';
else
    particleDiameter_m = table3ParticleSauterDiameter_m(params);
    if isfinite(particleDiameter_m) && particleDiameter_m > 0.0
        bedModel.particleDiameterSource = 'table3_d32_default';
    end
end
if ~(isfinite(particleDiameter_m) && particleDiameter_m > 0.0)
    particleDiameter_m = 0.015;
    bedModel.particleDiameterSource = 'fallback_constant';
end
bedModel.particleDiameter_m = particleDiameter_m;

fillHeight_m = max(params.bin.height_m * params.bin.fillFraction, 1e-12);
if isfield(params.aeration, 'representativeTubeCenterHeight_m') && ...
        ~isempty(params.aeration.representativeTubeCenterHeight_m) && ...
        isfinite(params.aeration.representativeTubeCenterHeight_m)
    repCenterHeight_m = params.aeration.representativeTubeCenterHeight_m;
else
    repCenterHeight_m = defaultRepresentativeTubeCenterHeight_m(fillHeight_m, tubeOD_m);
end
tubeRadius_m = max(0.5 * max(tubeOD_m, 0.0), 0.0);
repCenterHeight_m = min(max(repCenterHeight_m, tubeRadius_m), ...
    max(fillHeight_m - tubeRadius_m, tubeRadius_m));
verticalEscape_m = max(fillHeight_m - repCenterHeight_m, 1e-6);
dx_m = max(length_m, 1e-12) / max(N, 1);
nTubes = max(1, round(params.aeration.nParallelTubes));
tributaryArea_m2 = max(dx_m * max(params.bin.width_m, 1e-12) / max(nTubes, 1), 1e-12);
segmentCenters_m = dx_m * (((1:max(N, 1))' - 0.5));

bedModel.representativeTubeCenterHeight_m = repCenterHeight_m;
bedModel.verticalEscapeDistance_m = verticalEscape_m;
bedModel.tributaryArea_m2 = tributaryArea_m2;
bedModel.segmentCenters_m = segmentCenters_m;

openTop = isfield(params.bin, 'openTop') && logical(params.bin.openTop);
ventCount = 0;
ventDiameter_m = 0.0;
if isfield(params.bin, 'ventHoleCount') && ~isempty(params.bin.ventHoleCount)
    ventCount = max(0, round(params.bin.ventHoleCount));
end
if isfield(params.bin, 'ventHoleDiameter_m') && ~isempty(params.bin.ventHoleDiameter_m)
    ventDiameter_m = max(params.bin.ventHoleDiameter_m, 0.0);
end

if openTop
    bedModel.mode = 'open_top';
    bedModel.segmentPathLengths_m = verticalEscape_m * ones(max(N, 1), 1);
    bedModel.enabled = bedModel.porosity > 0 && bedModel.porosity < 1 && ...
        bedModel.particleDiameter_m > 0;
elseif ventCount > 0 && ventDiameter_m > 0
    bedModel.mode = 'covered_top_plus_vent';
    ventCenters_m = max(length_m, 1e-12) * (((1:ventCount)' - 0.5) / ventCount);
    axialOffsets_m = min(abs(segmentCenters_m - ventCenters_m'), [], 2);
    bedModel.segmentPathLengths_m = hypot(verticalEscape_m, axialOffsets_m);
    bedModel.enabled = bedModel.porosity > 0 && bedModel.porosity < 1 && ...
        bedModel.particleDiameter_m > 0;
    bedModel.derivedBoundaryEnabled = true;
    bedModel.outletAreaPerBranch_m2 = max( ...
        ventCount * pi * ventDiameter_m^2 / 4 / max(nTubes, 1), ...
        1e-12);
    bedModel.outletDischargeCoefficient = 0.63;
end
end

function boundaryGauge_Pa = derivedOutletBoundaryGaugePressure_Pa( ...
        branchMdot_kg_s, rhoBoundary_kg_m3, bedModel)
boundaryGauge_Pa = 0.0;
if isfield(bedModel, 'boundaryGaugeOffset_Pa') && ~isempty(bedModel.boundaryGaugeOffset_Pa)
    boundaryGauge_Pa = bedModel.boundaryGaugeOffset_Pa;
elseif isfield(bedModel, 'boundaryGaugePressure_Pa') && ~isempty(bedModel.boundaryGaugePressure_Pa)
    boundaryGauge_Pa = bedModel.boundaryGaugePressure_Pa;
end

if branchMdot_kg_s <= 0 || ~(isfield(bedModel, 'derivedBoundaryEnabled') && bedModel.derivedBoundaryEnabled)
    return;
end

area_m2 = NaN;
Cd = NaN;
if isfield(bedModel, 'outletAreaPerBranch_m2')
    area_m2 = bedModel.outletAreaPerBranch_m2;
end
if isfield(bedModel, 'outletDischargeCoefficient')
    Cd = bedModel.outletDischargeCoefficient;
end
if ~(isfinite(area_m2) && area_m2 > 0 && isfinite(Cd) && Cd > 0)
    return;
end

boundaryGauge_Pa = boundaryGauge_Pa + ...
    0.5 * (max(branchMdot_kg_s, 0.0) / max(Cd * area_m2, 1e-12))^2 / ...
    max(rhoBoundary_kg_m3, 1e-12);
end

function dp_Pa = ergunSegmentBedPressureDrop_Pa(mdotRelease_kg_s, rho_kg_m3, mu_Pa_s, bedModel, segmentIndex)
dp_Pa = 0.0;
if ~isfield(bedModel, 'enabled') || ~bedModel.enabled || mdotRelease_kg_s <= 0
    return;
end

porosity = bedModel.porosity;
eps3 = max(porosity^3, 1e-12);
oneMinusEps = max(1.0 - porosity, 0.0);
particleDiameter_m = max(bedModel.particleDiameter_m, 1e-12);
flowArea_m2 = max(bedModel.tributaryArea_m2, 1e-12);
idx = min(max(round(segmentIndex), 1), numel(bedModel.segmentPathLengths_m));
pathLength_m = max(bedModel.segmentPathLengths_m(idx), 0.0);
uSuperficial_m_s = max(mdotRelease_kg_s, 0.0) / ...
    max(max(rho_kg_m3, 1e-12) * flowArea_m2, 1e-12);
dpdx_Pa_m = ...
    150.0 * (oneMinusEps^2) / eps3 * max(mu_Pa_s, 0.0) * uSuperficial_m_s / (particleDiameter_m^2) + ...
    1.75 * oneMinusEps / eps3 * max(rho_kg_m3, 0.0) * uSuperficial_m_s^2 / particleDiameter_m;
dp_Pa = dpdx_Pa_m * pathLength_m;
end

function [mdotRelease_kg_s, bedGauge_Pa] = solveSegmentReleaseWithBedBackpressure( ...
        params, gaugePressure_Pa, rho_kg_m3, mu_Pa_s, totalHoleArea_m2, ...
        remainingMdot_kg_s, bedModel, segmentIndex, boundaryGauge_Pa)
mdotRelease_kg_s = 0.0;
bedGauge_Pa = boundaryGauge_Pa;
Cd = perforationDischargeCoefficient(params);

if remainingMdot_kg_s <= 0 || Cd <= 0 || totalHoleArea_m2 <= 0 || gaugePressure_Pa <= boundaryGauge_Pa
    return;
end

if ~(isfield(bedModel, 'enabled') && bedModel.enabled)
    mdotRelease_kg_s = min(perforationDischargeMassFlow_kg_s( ...
        params, gaugePressure_Pa, rho_kg_m3, totalHoleArea_m2), ...
        remainingMdot_kg_s);
    return;
end

    function [residual_kg_s, localBedGauge_Pa] = dischargeResidual(mdotTrial_kg_s)
        localBedGauge_Pa = boundaryGauge_Pa + ergunSegmentBedPressureDrop_Pa( ...
            mdotTrial_kg_s, rho_kg_m3, mu_Pa_s, bedModel, segmentIndex);
        predicted_kg_s = Cd * max(totalHoleArea_m2, 0.0) * ...
            sqrt(2 * max(rho_kg_m3, 0.0) * max(gaugePressure_Pa - localBedGauge_Pa, 0.0));
        residual_kg_s = predicted_kg_s - mdotTrial_kg_s;
    end

[residualHi_kg_s, bedGaugeHi_Pa] = dischargeResidual(remainingMdot_kg_s);
if residualHi_kg_s >= 0
    mdotRelease_kg_s = remainingMdot_kg_s;
    bedGauge_Pa = bedGaugeHi_Pa;
    return;
end

lo_kg_s = 0.0;
hi_kg_s = remainingMdot_kg_s;
for iter = 1:24
    mid_kg_s = 0.5 * (lo_kg_s + hi_kg_s);
    [residualMid_kg_s, bedGaugeMid_Pa] = dischargeResidual(mid_kg_s); %#ok<NASGU>
    if residualMid_kg_s >= 0
        lo_kg_s = mid_kg_s;
    else
        hi_kg_s = mid_kg_s;
    end
end
mdotRelease_kg_s = 0.5 * (lo_kg_s + hi_kg_s);
bedGauge_Pa = boundaryGauge_Pa + ergunSegmentBedPressureDrop_Pa( ...
    mdotRelease_kg_s, rho_kg_m3, mu_Pa_s, bedModel, segmentIndex);
end

function mdot_kg_s = perforationDischargeMassFlow_kg_s(params, gaugePressure_Pa, rho_kg_m3, totalHoleArea_m2)
deltaP_Pa = max(gaugePressure_Pa - aerationBedSideGaugePressure(params), 0.0);
mdot_kg_s = perforationDischargeCoefficient(params) * max(totalHoleArea_m2, 0.0) * ...
    sqrt(2 * max(rho_kg_m3, 0.0) * deltaP_Pa);
end

function [mEvap_kg_s, qEvap_W] = evaporationLossSegment( ...
        params, dmRelease_kg_s, Tair_C, Tbed_C, qAvail_W, holeVelocity_m_s)
surfaceT_C = min(max(params.evaporation.surfaceTemp_C, Tbed_C), Tair_C);
surfaceT_K = surfaceT_C + 273.15;
airInT_K = params.environment.greenhouseAir_C + 273.15;
props = airProps(max(Tair_C, params.environment.greenhouseAir_C) + 273.15, ...
    params.environment.pressure_Pa);

rhoV_surface_kg_m3 = waterVaporDensityAtSaturation(surfaceT_K);
rhoV_inlet_kg_m3 = params.evaporation.relativeHumidity * ...
    waterVaporDensityAtSaturation(airInT_K);
deltaRhoV_kg_m3 = max(0, rhoV_surface_kg_m3 - rhoV_inlet_kg_m3);

ReHole = props.rho_kg_m3 * max(holeVelocity_m_s, 0) * ...
    max(params.perforation.holeDiameter_m, 1e-9) / max(props.mu_Pa_s, 1e-12);
NuHole = churchillBernsteinNu(max(ReHole, 1e-9), props.Pr);
hEvap_W_m2K = params.calibration.evaporationHMultiplier * ...
    NuHole * props.k_W_mK / max(params.perforation.holeDiameter_m, 1e-9);

% [R10] Lewis relation for air-water systems with Lewis factor near unity:
% h_m = h / (rho*cp*Le^(2/3))
hMass_m_s = hEvap_W_m2K / max( ...
    props.rho_kg_m3 * props.cp_J_kgK * params.evaporation.lewisFactor^(2/3), ...
    1e-12);

segArea_m2 = params.evaporation.wettedArea_m2 / max( ...
    params.aeration.nParallelTubes * params.aeration.segments, 1);
volRelease_m3_s = dmRelease_kg_s / max(props.rho_kg_m3, 1e-12);
hfg_J_kg = latentHeatVaporizationWater(surfaceT_K);

% [R5,R10] Diffusion-limited evaporation from a wetted surface:
% m_dot = h_m * A * (rho_v,s - rho_v,inf)
mDiff_kg_s = hMass_m_s * segArea_m2 * deltaRhoV_kg_m3;

% Volumetric carrying-capacity upper bound for the released air parcel.
mCapacity_kg_s = volRelease_m3_s * deltaRhoV_kg_m3;

% [R9] Heat-limited evaporation using latent heat from the released-air
% sensible transfer available at this segment.
mHeat_kg_s = qAvail_W / max(hfg_J_kg, 1e-9);

mEvap_kg_s = max(0, min([mDiff_kg_s, mCapacity_kg_s, mHeat_kg_s]));
qEvap_W = mEvap_kg_s * hfg_J_kg;
end

function rhoV_kg_m3 = waterVaporDensityAtSaturation(T_K)
psat_Pa = saturationPressureWater_Pa(T_K);
Rv_J_kgK = 461.526;
rhoV_kg_m3 = psat_Pa / max(Rv_J_kgK * T_K, 1e-12);
end

function hfg_J_kg = latentHeatVaporizationWater(T_K)
rhoL_kg_m3 = saturationLiquidDensityWater(T_K);
rhoV_kg_m3 = saturationVaporDensityWater(T_K);
dpdT_Pa_K = saturationPressureSlopeWater(T_K);

% [R9] Clausius-Clapeyron relation on the saturation line:
% h_fg = T * (dp/dT) * (1/rho_v - 1/rho_l)
hfg_J_kg = T_K * dpdT_Pa_K * ...
    (1 / max(rhoV_kg_m3, 1e-12) - 1 / max(rhoL_kg_m3, 1e-12));
end

function psat_Pa = saturationPressureWater_Pa(T_K)
Tc_K = 647.096;
pc_Pa = 22.064e6;
theta = T_K / Tc_K;
tau = 1 - theta;
a = [-7.85951783 1.84408259 -11.7866497 22.6807411 -15.9618719 1.80122502];
exponents = [1 1.5 3 3.5 4 7.5];
series = 0;
for i = 1:numel(a)
    series = series + a(i) * tau^exponents(i);
end
% [R9] IAPWS saturation-pressure correlation, Eq. (1).
psat_Pa = pc_Pa * exp((Tc_K / T_K) * series);
end

function dpdT_Pa_K = saturationPressureSlopeWater(T_K)
Tc_K = 647.096;
pc_Pa = 22.064e6;
tau = 1 - T_K / Tc_K;
a = [-7.85951783 1.84408259 -11.7866497 22.6807411 -15.9618719 1.80122502];
exponents = [1 1.5 3 3.5 4 7.5];
sumSeries = 0;
sumDeriv = 0;
for i = 1:numel(a)
    sumSeries = sumSeries + a(i) * tau^exponents(i);
    sumDeriv = sumDeriv + a(i) * exponents(i) * tau^(exponents(i) - 1);
end
lnp = (Tc_K / T_K) * sumSeries;
psat_Pa = pc_Pa * exp(lnp);
dlndp_dT = -(Tc_K / T_K^2) * sumSeries - sumDeriv / T_K;
dpdT_Pa_K = psat_Pa * dlndp_dT;
end

function rhoL_kg_m3 = saturationLiquidDensityWater(T_K)
Tc_K = 647.096;
rhoc_kg_m3 = 322;
tau = 1 - T_K / Tc_K;
b = [1.99274064 1.09965342 -0.510839303 -1.75493479 -45.5170352 -6.74694450e5];
exponents = [1/3 2/3 5/3 16/3 43/3 110/3];
series = 1;
for i = 1:numel(b)
    series = series + b(i) * tau^exponents(i);
end
% [R9] IAPWS saturated-liquid density correlation, Eq. (2).
rhoL_kg_m3 = rhoc_kg_m3 * series;
end

function rhoV_kg_m3 = saturationVaporDensityWater(T_K)
Tc_K = 647.096;
rhoc_kg_m3 = 322;
tau = 1 - T_K / Tc_K;
c = [-2.03150240 -2.68302940 -5.38626492 -17.2991605 -44.7586581 -63.9201063];
exponents = [2/6 4/6 8/6 18/6 37/6 71/6];
series = 0;
for i = 1:numel(c)
    series = series + c(i) * tau^exponents(i);
end
% [R9] IAPWS saturated-vapor density correlation, Eq. (3).
rhoV_kg_m3 = rhoc_kg_m3 * exp(series);
end

function dp_Pa = abruptExpansionLoss(mdot_kg_s, T_K, D1_m, D2_m, pressure_Pa)
props = airProps(T_K, pressure_Pa);
A1_m2 = pi * D1_m^2 / 4;
A2_m2 = pi * D2_m^2 / 4;
u1_m_s = mdot_kg_s / max(props.rho_kg_m3 * A1_m2, 1e-12);
epsilon = A1_m2 / max(A2_m2, 1e-12);

% [R6] Idel'chik Section IV sudden expansion; identical to the standard
% Borda-Carnot coefficient on the upstream-velocity basis:
% zeta = (1 - A1/A2)^2.
K = max((1 - epsilon)^2, 0);
dp_Pa = K * 0.5 * props.rho_kg_m3 * u1_m_s^2;
end

function dp_Pa = suddenContractionLoss(mdot_kg_s, T_K, D1_m, D2_m, pressure_Pa, angle_deg)
if D2_m >= D1_m
    dp_Pa = 0;
    return;
end
props = airProps(T_K, pressure_Pa);
A1_m2 = pi * D1_m^2 / 4;
A2_m2 = pi * D2_m^2 / 4;
epsilon = max(min(A2_m2 / max(A1_m2, 1e-12), 1), 1e-12);
u2_m_s = mdot_kg_s / max(props.rho_kg_m3 * A2_m2, 1e-12);

% [R12] Idel'chik, Section III, Diagram 3-5: conical converging bellmouth
% without end wall. The tabulated coefficient C is used with Eq. (3-3):
% zeta = C * (1/epsilon - 1), referenced to the smaller-pipe velocity.
C = idelchikConicalBellmouthC(epsilon, angle_deg);
K = max(C * (1 / max(epsilon, 1e-12) - 1), 0);
dp_Pa = K * 0.5 * props.rho_kg_m3 * u2_m_s^2;
end

function C = idelchikConicalBellmouthC(epsilon, angle_deg)
epsTable = [0.025 0.050 0.075 0.10 0.15 0.25 0.60 1.0];
angleTable_deg = [0 10 20 30 40 60 100 140 180];
CTable = [ ...
    1.00 0.96 0.93 0.90 0.86 0.80 0.69 0.59 0.50; ...
    1.00 0.93 0.86 0.80 0.75 0.67 0.58 0.53 0.50; ...
    1.00 0.87 0.75 0.65 0.58 0.50 0.48 0.49 0.50; ...
    1.00 0.80 0.69 0.55 0.48 0.41 0.41 0.44 0.50; ...
    1.00 0.76 0.58 0.43 0.33 0.25 0.27 0.38 0.50; ...
    1.00 0.68 0.45 0.30 0.22 0.17 0.22 0.34 0.50; ...
    1.00 0.46 0.27 0.18 0.14 0.13 0.21 0.33 0.50; ...
    1.00 0.32 0.20 0.14 0.11 0.10 0.18 0.30 0.50];

epsClamped = max(min(epsilon, epsTable(end)), epsTable(1));
angleClamped = max(min(angle_deg, angleTable_deg(end)), angleTable_deg(1));
colInterp = zeros(size(epsTable));
for i = 1:numel(epsTable)
    colInterp(i) = interp1(angleTable_deg, CTable(i, :), angleClamped, 'linear');
end
C = interp1(epsTable, colInterp, epsClamped, 'linear');
end

function K = hooper2KFittingLossCoefficient(Re, diameter_m, K1, Kinf)
if Re <= 1e-12 || diameter_m <= 0
    K = 0;
    return;
end
diameter_in = max(diameter_m * 39.37007874015748, 1e-12);
% [R14] Crowl and Louvar, Eq. (4-38): Hooper 2-K fitting-loss method.
K = K1 / Re + Kinf * (1 + 1 / diameter_in);
end

function K = hooper2KEntranceLossCoefficient(Re, K1, Kinf)
if nargin < 2
    K1 = 160;
end
if nargin < 3
    Kinf = 0.50;
end
if Re <= 1e-12
    K = 0;
    return;
end
% [R14] Crowl and Louvar, Eq. (4-39): entrance/exit form of the 2-K method.
K = K1 / Re + Kinf;
end

function multiplier = aerationLossMultiplier(params, fieldName, legacyFieldName)
if isfield(params.aeration, fieldName)
    multiplier = max(0, params.aeration.(fieldName));
elseif isfield(params.aeration, legacyFieldName)
    multiplier = max(0, params.aeration.(legacyFieldName));
else
    multiplier = 1.0;
end
end

function losses = distributionNetworkLosses(params, totalFlow_Lpm, branchFlow_Lpm, T_K, downstreamID_m)
props = airProps(T_K, params.environment.pressure_Pa);
volFlowTotal_m3_s = totalFlow_Lpm / (1000 * 60);
volFlowBranch_m3_s = branchFlow_Lpm / (1000 * 60);
mdotBranch_kg_s = props.rho_kg_m3 * volFlowBranch_m3_s;
splitterOutletCount = max(1, round(params.aeration.splitterOutletCount));
activeSplitterBranches = min(splitterOutletCount, max(1, params.aeration.nParallelTubes));
mdotSplitter_kg_s = mdotBranch_kg_s * activeSplitterBranches;

connectorID_m = params.aeration.branchConnectorID_m;
connectorArea_m2 = pi * connectorID_m^2 / 4;
uConnector_m_s = mdotBranch_kg_s / max(props.rho_kg_m3 * connectorArea_m2, 1e-12);
ReConnector = props.rho_kg_m3 * uConnector_m_s * connectorID_m / ...
    max(props.mu_Pa_s, 1e-12);
fConnector = churchillFrictionFactor(ReConnector);
headerLossMultiplier = aerationLossMultiplier(params, 'headerLossMultiplier', 'headerMinorK');
splitterBodyLossMultiplier = aerationLossMultiplier(params, 'splitterBodyLossMultiplier', 'splitterBodyK');
splitterBranchTeeK = hooper2KFittingLossCoefficient(ReConnector, connectorID_m, 500, 0.70);
splitterValveK = hooper2KFittingLossCoefficient(ReConnector, connectorID_m, 300, 0.10);

losses = struct();
losses.header_Pa = 0;
losses.headerToSplitter_Pa = 0;
losses.headerLossCoefficient = 0;
losses.headerLossMultiplier = headerLossMultiplier;
losses.splitterBranchTeeLossCoefficient = splitterBranchTeeK;
losses.splitterValveLossCoefficient = splitterValveK;
losses.splitterBodyLossCoefficient = splitterBodyLossMultiplier * ...
    (splitterBranchTeeK + splitterValveK);
losses.splitterBodyLossMultiplier = splitterBodyLossMultiplier;
losses.splitterBody_Pa = losses.splitterBodyLossCoefficient * ...
    0.5 * props.rho_kg_m3 * uConnector_m_s^2;
losses.connectorFriction_Pa = fConnector * ...
    params.aeration.branchConnectorLength_m / max(connectorID_m, 1e-12) * ...
    0.5 * props.rho_kg_m3 * uConnector_m_s^2;
losses.connectorToBranch_Pa = 0;

if isfield(params.aeration, 'headerEnabled') && params.aeration.headerEnabled
    Aheader_m2 = pi * params.aeration.headerID_m^2 / 4;
    uHeader_m_s = volFlowTotal_m3_s / max(Aheader_m2, 1e-12);
    ReHeader = props.rho_kg_m3 * uHeader_m_s * params.aeration.headerID_m / ...
        max(props.mu_Pa_s, 1e-12);

    % [R14] The unresolved optional upstream-header entry is represented
    % by the normal-entrance 2-K correlation, optionally scaled by a
    % calibration multiplier.
    losses.headerLossCoefficient = headerLossMultiplier * ...
        hooper2KEntranceLossCoefficient(ReHeader, 160, 0.50);
    losses.header_Pa = losses.headerLossCoefficient * ...
        0.5 * props.rho_kg_m3 * uHeader_m_s^2;
    losses.headerToSplitter_Pa = suddenContractionLoss( ...
        mdotSplitter_kg_s, T_K, params.aeration.headerID_m, ...
        params.aeration.splitterInletID_m, params.environment.pressure_Pa, ...
        params.aeration.contractionConeAngle_deg);
end

if downstreamID_m < connectorID_m
    losses.connectorToBranch_Pa = suddenContractionLoss( ...
        mdotBranch_kg_s, T_K, connectorID_m, downstreamID_m, ...
        params.environment.pressure_Pa, ...
        params.aeration.contractionConeAngle_deg);
elseif downstreamID_m > connectorID_m
    losses.connectorToBranch_Pa = abruptExpansionLoss( ...
        mdotBranch_kg_s, T_K, connectorID_m, downstreamID_m, ...
        params.environment.pressure_Pa);
end

losses.total_Pa = losses.header_Pa + losses.headerToSplitter_Pa + ...
    losses.splitterBody_Pa + losses.connectorFriction_Pa + ...
    losses.connectorToBranch_Pa;
end

function [Tb_C, Tt_C] = twoNodeSteadyState(binModel, Tamb_C, Q_W, bottomFraction)
Uab = binModel.UAbottomAmbient_W_K;
Uat = binModel.UAtopAmbient_W_K;
U12 = binModel.UAinternal_W_K;
Qb = bottomFraction * Q_W;
Qt = (1 - bottomFraction) * Q_W;

% [R5] Steady energy balance on the two bed nodes.
A = [Uab + U12, -U12; -U12, Uat + U12];
b = [Qb + Uab * Tamb_C; Qt + Uat * Tamb_C];
x = A \ b;

Tb_C = x(1);
Tt_C = x(2);
end

function transient = simulateBinPID(params, binModel, recommended)
dt_s = params.control.dt_s;
duration_s = params.control.duration_h * 3600;
t_s = (0:dt_s:duration_s).';
n = numel(t_s);

Tb_C = zeros(n, 1);
Tt_C = zeros(n, 1);
Tamb_C = params.environment.greenhouseAir_C * ones(n, 1);
Vcmd_V = zeros(n, 1);
Qheater_W = zeros(n, 1);
TairOut_C = zeros(n, 1);
TwireMax_C = zeros(n, 1);

Tb_C(1) = params.control.initialBottom_C;
Tt_C(1) = params.control.initialTop_C;

integralError = 0;
prevError = params.bin.setpoint_C - min(Tb_C(1), Tt_C(1));

for k = 2:n
    Tcontrol_C = min(Tb_C(k - 1), Tt_C(k - 1));
    error_C = params.bin.setpoint_C - Tcontrol_C;
    derivative_C_s = (error_C - prevError) / dt_s;

    tentativeIntegral = integralError + error_C * dt_s;
    unsatVoltage_V = params.control.Kp * error_C + ...
        params.control.Ki * tentativeIntegral + ...
        params.control.Kd * derivative_C_s;

    Vmax_V = recommended.voltage_V;
    Vcmd_V(k) = min(max(unsatVoltage_V, 0), Vmax_V);

    if (Vcmd_V(k) > 0 && Vcmd_V(k) < Vmax_V) || ...
            (Vcmd_V(k) == 0 && error_C > 0) || ...
            (Vcmd_V(k) == Vmax_V && error_C < 0)
        integralError = tentativeIntegral;
    end

    if max(Tb_C(k - 1), Tt_C(k - 1)) >= params.control.highCutout_C
        Vcmd_V(k) = 0;
    end

    Tmean_C = 0.5 * (Tb_C(k - 1) + Tt_C(k - 1));
    op = evaluateOperatingPoint( ...
        params, Vcmd_V(k), recommended.totalFlow_Lpm, Tmean_C, Tamb_C(k));

    Qheater_W(k) = op.total.QtoBed_W;
    TairOut_C(k) = op.branch.heater.airOutletTemp_C;
    TwireMax_C(k) = op.branch.heater.wireMaxTemp_C;

    Qbottom_W = params.bin.heatToBottomFraction * Qheater_W(k);
    Qtop_W = (1 - params.bin.heatToBottomFraction) * Qheater_W(k);

    % [R5] Transient two-node energy balances.
    dTb_dt_K_s = ( ...
        Qbottom_W + ...
        binModel.UAinternal_W_K * (Tt_C(k - 1) - Tb_C(k - 1)) - ...
        binModel.UAbottomAmbient_W_K * (Tb_C(k - 1) - Tamb_C(k))) / ...
        max(binModel.Cbottom_J_K, 1e-9);

    dTt_dt_K_s = ( ...
        Qtop_W + ...
        binModel.UAinternal_W_K * (Tb_C(k - 1) - Tt_C(k - 1)) - ...
        binModel.UAtopAmbient_W_K * (Tt_C(k - 1) - Tamb_C(k))) / ...
        max(binModel.Ctop_J_K, 1e-9);

    Tb_C(k) = Tb_C(k - 1) + dTb_dt_K_s * dt_s;
    Tt_C(k) = Tt_C(k - 1) + dTt_dt_K_s * dt_s;
    prevError = error_C;
end

transient = struct();
transient.time_h = t_s / 3600;
transient.ambient_C = Tamb_C;
transient.bottom_C = Tb_C;
transient.top_C = Tt_C;
transient.voltage_V = Vcmd_V;
transient.heatToBed_W = Qheater_W;
transient.airOutlet_C = TairOut_C;
transient.wireMax_C = TwireMax_C;
transient.onFraction = mean(Vcmd_V > 0.1);
end

function props = airProps(T_K, pressure_Pa)
T_K = max(T_K, 200);
R_J_kgK = 287.058;
cp_J_kgK = 1007;
muRef_Pa_s = 1.716e-5;
Tref_K = 273.15;
S_K = 111.0;

% [R7] Sutherland law for air viscosity.
mu_Pa_s = muRef_Pa_s * (T_K / Tref_K)^(3/2) * (Tref_K + S_K) / (T_K + S_K);

Pr = 0.71;

% [R7] k = mu*cp/Pr with constant Pr approximation for air.
k_W_mK = mu_Pa_s * cp_J_kgK / Pr;

% [R7] rho = p / (R*T)
rho_kg_m3 = pressure_Pa / (R_J_kgK * T_K);

props = struct();
props.cp_J_kgK = cp_J_kgK;
props.mu_Pa_s = mu_Pa_s;
props.k_W_mK = k_W_mK;
props.rho_kg_m3 = rho_kg_m3;
props.Pr = Pr;
end

function wireProps = nichrome80Props(wire, T_C)
area_m2 = pi * wire.diameter_m^2 / 4;

% Linearized temperature dependence for resistivity near the design range.
rho_Ohm_m = wire.rho20_Ohm_m * (1 + wire.tempCoeff_1_K * (T_C - 20));

wireProps = struct();
wireProps.area_m2 = area_m2;
wireProps.rho_Ohm_m = rho_Ohm_m;
wireProps.resistance_Ohm_m = rho_Ohm_m / max(area_m2, 1e-12);
end

function Nu = churchillBernsteinNu(Re, Pr)
Re = max(Re, 1e-9);
Pr = max(Pr, 1e-9);

% [R1] Churchill-Bernstein crossflow correlation for a cylinder.
Nu = 0.3 + ...
    (0.62 * Re^0.5 * Pr^(1/3)) / (1 + (0.4 / Pr)^(2/3))^(1/4) * ...
    (1 + (Re / 282000)^(5/8))^(4/5);
end

function Nu = internalTubeNu(Re, Pr, D_m, L_m)
Re = max(Re, 1e-9);
Pr = max(Pr, 1e-9);

if Re < 2300
    % [R5] Hausen-type developing-laminar relation via Graetz number.
    Gz = Re * Pr * D_m / max(L_m, 1e-9);
    Nu = 3.66 + (0.0668 * Gz) / (1 + 0.04 * Gz^(2/3));
else
    % [R2] Gnielinski turbulent internal-flow correlation.
    f = churchillFrictionFactor(Re);
    Nu = (f / 8) * (Re - 1000) * Pr / ...
        (1 + 12.7 * sqrt(f / 8) * (Pr^(2/3) - 1));
    Nu = max(Nu, 3.66);
end
end

function f = churchillFrictionFactor(Re)
Re = max(Re, 1e-9);
if Re < 2300
    % [R5] Laminar pipe flow: f = 64 / Re
    f = 64 / Re;
else
    % [R3] Churchill explicit smooth-tube friction-factor form.
    A = (2.457 * log(1 / ((7 / Re)^0.9)))^16;
    B = (37530 / Re)^16;
    f = 8 * ((8 / Re)^12 + 1 / (A + B)^(1.5))^(1 / 12);
end
end

function printSummary(results)
 [activeCfg, activeScenario, rec] = comparisonBasedReference(results);

fprintf('\n============================================================\n');
fprintf('VERMICOMPOSTER AIR-HEATER DESIGN MODEL\n');
fprintf('============================================================\n');
fprintf('Greenhouse air design temperature          : %.1f C\n', ...
    activeScenario.params.environment.greenhouseAir_C);
fprintf('Active vessel configuration                : %s\n', ...
    activeScenario.binModel.topMode);
if ~activeScenario.params.bin.openTop && activeScenario.binModel.ventHoleArea_m2 > 0
    fprintf('Top vent hole(s)                           : %d x %.1f mm dia (total open area %.6f m^2)\n', ...
        activeScenario.params.bin.ventHoleCount, ...
        1000 * activeScenario.params.bin.ventHoleDiameter_m, ...
        activeScenario.binModel.ventHoleArea_m2);
end
if ~isempty(fieldnames(activeCfg))
    fprintf('Comparison-plot reference row              : %s\n', ...
        activeCfg.label);
end
fprintf('\nDesign specifications:\n');
fprintf('  Vessel interior                         : %.2f m x %.2f m x %.2f m, fill fraction %.2f\n', ...
    activeScenario.params.bin.length_m, activeScenario.params.bin.width_m, ...
    activeScenario.params.bin.height_m, activeScenario.params.bin.fillFraction);
fprintf('  Vessel wall stack                       : %.1f mm steel (k = %.1f W/m-K) + %.1f mm insulation (k = %.3f W/m-K)\n', ...
    1000 * activeScenario.params.binWall.sheetThickness_m, ...
    activeScenario.params.binWall.sheetK_W_mK, ...
    1000 * activeScenario.params.binWall.insulationThickness_m, ...
    activeScenario.params.binWall.insulationK_W_mK);
fprintf('  Heater tube                             : L %.3f m, ID %.1f mm, OD %.1f mm, insulation %.1f mm\n', ...
    activeScenario.params.heaterTube.length_m, ...
    1000 * activeScenario.params.heaterTube.ID_m, ...
    1000 * activeScenario.params.heaterTube.OD_m, ...
    1000 * activeScenario.params.heaterTube.insulationThickness_m);
fprintf('  Heater coil geometry                    : mean D %.1f mm, pitch %.1f mm, axial span %.1f mm\n', ...
    1000 * activeScenario.params.heaterTube.coilMeanD_m, ...
    rec.pitch_mm, ...
    1000 * activeScenario.params.heaterTube.coilSpan_m);
fprintf('  Aeration network                        : %d tubes, each L %.2f m, ID %.1f mm, OD %.1f mm\n', ...
    activeScenario.params.aeration.nParallelTubes, ...
    activeScenario.params.aeration.length_m, ...
    1000 * activeScenario.params.aeration.ID_m, ...
    1000 * activeScenario.params.aeration.OD_m);
if isfield(activeScenario.params.aeration, 'headerEnabled') && activeScenario.params.aeration.headerEnabled
    fprintf('  Distribution manifold                  : header ID %.1f mm feeding %d-way splitter(s), splitter inlet ID %.1f mm, branch connector ID %.1f mm, connector length %.3f m\n', ...
        1000 * activeScenario.params.aeration.headerID_m, ...
        activeScenario.params.aeration.splitterOutletCount, ...
        1000 * activeScenario.params.aeration.splitterInletID_m, ...
        1000 * activeScenario.params.aeration.branchConnectorID_m, ...
        activeScenario.params.aeration.branchConnectorLength_m);
else
    fprintf('  Distribution manifold                  : no separate header; %d-way splitter(s), splitter inlet ID %.1f mm, branch connector ID %.1f mm, connector length %.3f m\n', ...
        activeScenario.params.aeration.splitterOutletCount, ...
        1000 * activeScenario.params.aeration.splitterInletID_m, ...
        1000 * activeScenario.params.aeration.branchConnectorID_m, ...
        activeScenario.params.aeration.branchConnectorLength_m);
end
fprintf('  Contraction model                      : Idel''chik conical converging bellmouth without end wall, %.0f deg\n', ...
    activeScenario.params.aeration.contractionConeAngle_deg);
fprintf('  Heater wire                             : %s, %s, %.3f mm diameter\n', ...
    activeScenario.params.wire.name, activeScenario.params.wire.gauge, ...
    1000 * activeScenario.params.wire.diameter_m);
fprintf('  Operating sweep                         : V = %.0f:%.0f:%.0f V, flow = %.0f:%.0f:%.0f L/min\n', ...
    activeScenario.params.sweep.voltage_V(1), ...
    activeScenario.params.sweep.voltage_V(2) - activeScenario.params.sweep.voltage_V(1), ...
    activeScenario.params.sweep.voltage_V(end), ...
    activeScenario.params.sweep.totalFlow_Lpm(1), ...
    activeScenario.params.sweep.totalFlow_Lpm(2) - activeScenario.params.sweep.totalFlow_Lpm(1), ...
    activeScenario.params.sweep.totalFlow_Lpm(end));
fprintf('  Limits                                  : I_total <= %.1f A, T_wire <= %.1f C, T_air,out <= %.1f C, dP <= %.0f Pa, V_hole <= %.2f m/s\n', ...
    activeScenario.params.limits.maxTotalCurrent_A, ...
    activeScenario.params.limits.maxWireTemp_C, ...
    activeScenario.params.limits.maxAirOutletTemp_C, ...
    activeScenario.params.limits.maxPressureDrop_Pa, ...
    activeScenario.params.limits.maxHoleVelocity_m_s);
fprintf('Lumped heat-loss requirement at %.1f C bed  : %.1f W\n', ...
    activeScenario.params.bin.setpoint_C, activeScenario.requiredHeat_W);
fprintf('Overall UA to greenhouse air               : %.3f W/K\n', ...
    activeScenario.lumpedLoss.UA_W_K);
fprintf('External natural-convection h estimate     : %.2f W/m^2-K\n', ...
    activeScenario.binModel.externalH_W_m2K);
fprintf('Lumped time constant                       : %.1f h\n', ...
    activeScenario.lumpedLoss.tau_h);
fprintf('Lumped Biot number estimate                : %.2f\n', ...
    activeScenario.lumpedLoss.Bi_lumped);
if activeScenario.lumpedLoss.lumpedStrictlyValid
    fprintf('Lumped-capacitance criterion               : Bi < 0.1 satisfied\n');
else
    fprintf('Lumped-capacitance criterion               : Bi < 0.1 NOT satisfied\n');
    fprintf('                                             use the lumped model for\n');
    fprintf('                                             heat-loss sizing, not for\n');
    fprintf('                                             detailed internal gradients\n');
end

fprintf('Minimum total flow at 100%% duty, %.1f C cap : %.1f L/min\n', ...
    activeScenario.params.limits.maxAirOutletTemp_C, ...
    activeScenario.minimumFlowAt100Duty_Lpm);
fprintf('Minimum total flow at %.0f%% duty, %.1f C cap : %.1f L/min\n', ...
    100 * activeScenario.params.sweep.targetDuty, ...
    activeScenario.params.limits.maxAirOutletTemp_C, ...
    activeScenario.minimumFlowAtTargetDuty_Lpm);
fprintf('Maximum allowed total current              : %.1f A\n', ...
    activeScenario.params.limits.maxTotalCurrent_A);

numFeasible = sum([activeScenario.designPoints.isFeasible]);
fprintf('Feasible operating points found            : %d of %d\n', ...
    numFeasible, numel(activeScenario.designPoints));

fprintf('\nVessel configuration comparison:\n');
fprintf('  config                    UA(W/K)  Qreq(W)  Flow100  Flow20  Pref(W)  Qref(L/min)  Iref(A)  status\n');
for i = 1:numel(results.configurationComparisons)
    scenario = results.configurationComparisons(i).scenario;
    recCfg = scenarioReferencePoint(scenario);
    fprintf('  %-24s %7.3f  %7.1f  %7.1f %7.1f %8.1f %11.1f %8.2f  %s\n', ...
        results.configurationComparisons(i).label, ...
        scenario.lumpedLoss.UA_W_K, ...
        scenario.requiredHeat_W, ...
        scenario.minimumFlowAt100Duty_Lpm, ...
        scenario.minimumFlowAtTargetDuty_Lpm, ...
        recCfg.totalPower_W, ...
        recCfg.totalFlow_Lpm, ...
        recCfg.totalCurrent_A, ...
        recommendedStatusShort(recCfg, scenario.params));
end

fprintf('\nCriteria-based recommended design point:\n');
if ~isempty(fieldnames(activeCfg))
    fprintf('  Configuration                      : %s\n', activeCfg.label);
end
fprintf('  Pitch                              : %.1f mm\n', rec.pitch_mm);
fprintf('  Supply voltage per branch          : %.1f V\n', rec.voltage_V);
fprintf('  Total flow                         : %.1f L/min\n', rec.totalFlow_Lpm);
fprintf('  Selected operating trio            : %.1f W, %.2f A, %.1f L/min\n', ...
    rec.totalPower_W, rec.totalCurrent_A, rec.totalFlow_Lpm);
fprintf('  Branch flow                        : %.1f L/min\n', rec.branchFlow_Lpm);
fprintf('  Branch wire length                 : %.2f m\n', rec.wireLength_m);
fprintf('  Branch coil axial length           : %.2f m\n', ...
    activeScenario.params.heaterTube.coilSpan_m);
fprintf('  Wire gauge / diameter              : %s / %.3f mm\n', ...
    activeScenario.params.wire.gauge, 1000 * activeScenario.params.wire.diameter_m);
if isfield(rec, 'electricalTopologyLabel')
    fprintf('  Electrical topology                : %s\n', rec.electricalTopologyLabel);
end
if isfield(rec, 'stringTubeCounts') && ~isempty(rec.stringTubeCounts)
    fprintf('  Series tube counts per string      : [%s]\n', ...
        sprintf('%g ', rec.stringTubeCounts));
end
fprintf('  Tube resistance                    : %.2f ohm\n', rec.branchResistance_Ohm);
fprintf('  Supply equivalent resistance       : %.2f ohm\n', rec.supplyEquivalentResistance_Ohm);
if isfield(rec, 'tubeCurrentMin_A')
    fprintf('  Tube current min/mean/max          : %.2f / %.2f / %.2f A\n', ...
        rec.tubeCurrentMin_A, rec.tubeCurrentMean_A, rec.tubeCurrentMax_A);
else
    fprintf('  Branch current                     : %.2f A\n', rec.branchCurrent_A);
end
fprintf('  Total current                      : %.2f A\n', rec.totalCurrent_A);
if isfield(rec, 'tubePowerMin_W')
    fprintf('  Tube power min/mean/max            : %.1f / %.1f / %.1f W\n', ...
        rec.tubePowerMin_W, rec.tubePowerMean_W, rec.tubePowerMax_W);
else
    fprintf('  Branch electrical power            : %.1f W\n', rec.branchPower_W);
end
fprintf('  Total electrical power             : %.1f W\n', rec.totalPower_W);
fprintf('  h_wire                             : %.1f W/m^2-K\n', rec.hWire_W_m2K);
fprintf('  Branch air outlet temperature      : %.1f C\n', rec.airOutlet_C);
fprintf('  Max wire temperature               : %.1f C\n', rec.wireMax_C);
fprintf('  Total heat to bed                  : %.1f W\n', rec.QtoBed_W);
fprintf('  Duty needed at design loss         : %.2f\n', rec.dutyNeeded);
fprintf('  Total pressure drop                : %.1f Pa\n', rec.deltaP_Pa);
if isfield(rec, 'distributionDeltaP_Pa')
    fprintf('  Distribution manifold dP          : %.1f Pa\n', rec.distributionDeltaP_Pa);
    fprintf('  Distribution dP components         : header %.1f, header->splitter %.1f, splitter body %.1f, connector friction %.1f, connector->branch %.1f Pa\n', ...
        rec.distributionHeader_Pa, rec.distributionHeaderToSplitter_Pa, ...
        rec.distributionSplitterBody_Pa, rec.distributionConnectorFriction_Pa, ...
        rec.distributionConnectorToBranch_Pa);
end
fprintf('  Full-power bed estimate            : bottom %.1f C, top %.1f C\n', ...
    rec.TbottomFullPower_C, rec.TtopFullPower_C);
if rec.isFeasible
    fprintf('  Status                             : meets all explicit limits and the heating target\n');
elseif explicitlySafe(rec, activeScenario.params)
    fprintf('  Status                             : meets explicit limits but does not fully meet the heating target\n');
else
    fprintf('  Status                             : no candidate in sweep satisfies all explicit limits and the heating target\n');
end

if rec.dutyNeeded > activeScenario.params.sweep.targetDuty
    fprintf('  Note                               : %.0f%% duty is not achieved;\n', ...
        100 * activeScenario.params.sweep.targetDuty);
    fprintf('                                       more tubes help local temperature\n');
    fprintf('                                       distribution, but not the minimum\n');
    fprintf('                                       m_dot*cp*DeltaT requirement\n');
end

topList = topCandidates(activeScenario.designPoints, 8, activeScenario.params);
fprintf('\nTop candidate points for active configuration:\n');
fprintf('  pitch  V   Qtot(L/min)  Ptot(W)  Itot(A)  Tair(C)  Twire(C)  Qbed(W)  duty   dP(Pa)\n');
for i = 1:numel(topList)
    fprintf('  %5.1f %2.0f    %8.1f   %7.1f   %7.2f   %7.1f   %8.1f   %7.1f  %4.2f   %6.1f\n', ...
        topList(i).pitch_mm, topList(i).voltage_V, topList(i).totalFlow_Lpm, ...
        topList(i).totalPower_W, topList(i).totalCurrent_A, ...
        topList(i).airOutlet_C, topList(i).wireMax_C, ...
        topList(i).QtoBed_W, topList(i).dutyNeeded, topList(i).deltaP_Pa);
end

if ~isempty(fieldnames(activeScenario.transient))
    fprintf('\nPID transient with recommended point as Vmax:\n');
    fprintf('  Final bottom temp                  : %.2f C\n', activeScenario.transient.bottom_C(end));
    fprintf('  Final top temp                     : %.2f C\n', activeScenario.transient.top_C(end));
    fprintf('  Peak wire temp                     : %.2f C\n', max(activeScenario.transient.wireMax_C));
    fprintf('  Peak air outlet                    : %.2f C\n', max(activeScenario.transient.airOutlet_C));
    fprintf('  Simulated heater on-fraction       : %.2f\n', activeScenario.transient.onFraction);
end

fprintf('============================================================\n\n');
end

function topList = topCandidates(designPoints, nTop, params)
criteria = criteriaMatrix(designPoints, params);
[~, order] = sortrows(criteria);
topList = designPoints(order(1:min(nTop, numel(designPoints))));
end

function txt = recommendedStatusShort(rec, params)
if isempty(fieldnames(rec))
    txt = 'no point';
elseif rec.isFeasible
    txt = 'feasible';
elseif explicitlySafe(rec, params)
    txt = 'hard-safe';
else
    txt = 'best-available only';
end
end

function tf = explicitlySafe(rec, params)
spread_C = abs(rec.TbottomFullPower_C - rec.TtopFullPower_C);
tf = rec.totalCurrent_A <= params.limits.maxTotalCurrent_A && ...
    rec.wireMax_C <= params.limits.maxWireTemp_C && ...
    rec.airOutlet_C <= params.limits.maxAirOutletTemp_C && ...
    rec.deltaP_Pa <= params.limits.maxPressureDrop_Pa && ...
    rec.waterLoss_kg_day <= params.limits.maxWaterLoss_kg_day && ...
    rec.holeVelocity_m_s <= params.limits.maxHoleVelocity_m_s && ...
    rec.TbottomFullPower_C >= params.optimization.minBottomTemp_C && ...
    rec.TtopFullPower_C >= params.optimization.minTopTemp_C && ...
    spread_C <= params.optimization.spreadLimit_C && ...
    rec.QtoBed_W > 0;
end

function gauge = nearestAwgLabel(diameter_m)
awgVals = 10:40;
diameters_m = 0.127e-3 * 92.^((36 - awgVals) / 39);
[~, idx] = min(abs(diameters_m - diameter_m));
gauge = sprintf('%d AWG', awgVals(idx));
end

function stateCode = feasibilityStateCode(pt, params)
if pt.isFeasible
    stateCode = 1;
elseif pt.totalCurrent_A > params.limits.maxTotalCurrent_A
    stateCode = 2;
elseif pt.airOutlet_C > params.limits.maxAirOutletTemp_C
    stateCode = 3;
elseif pt.wireMax_C > params.limits.maxWireTemp_C
    stateCode = 4;
elseif pt.deltaP_Pa > params.limits.maxPressureDrop_Pa
    stateCode = 5;
elseif pt.waterLoss_kg_day > params.limits.maxWaterLoss_kg_day
    stateCode = 6;
elseif pt.holeVelocity_m_s > params.limits.maxHoleVelocity_m_s
    stateCode = 7;
else
    stateCode = 8;
end
end

function makePlots(results)
[activeCfg, activeScenario, rec] = comparisonBasedReference(results);

pts = activeScenario.designPoints;
flowVals = unique([pts.totalFlow_Lpm]);
voltVals = unique([pts.voltage_V]);

if ~isempty(fieldnames(rec))
    baselinePitch = rec.pitch_mm;
else
    baselinePitch = activeScenario.params.sweep.pitch_mm(1);
end
sel = [pts.pitch_mm] == baselinePitch;
ptsBase = pts(sel);

Qgrid = nan(numel(flowVals), numel(voltVals));
TairGrid = nan(numel(flowVals), numel(voltVals));
DutyGrid = nan(numel(flowVals), numel(voltVals));
CurrentGrid = nan(numel(flowVals), numel(voltVals));
WireGrid = nan(numel(flowVals), numel(voltVals));
DPgrid = nan(numel(flowVals), numel(voltVals));
StateGrid = nan(numel(flowVals), numel(voltVals));

for i = 1:numel(ptsBase)
    row = find(flowVals == ptsBase(i).totalFlow_Lpm, 1, 'first');
    col = find(voltVals == ptsBase(i).voltage_V, 1, 'first');
    Qgrid(row, col) = ptsBase(i).QtoBed_W;
    TairGrid(row, col) = ptsBase(i).airOutlet_C;
    DutyGrid(row, col) = ptsBase(i).dutyNeeded;
    CurrentGrid(row, col) = ptsBase(i).totalCurrent_A;
    WireGrid(row, col) = ptsBase(i).wireMax_C;
    DPgrid(row, col) = ptsBase(i).deltaP_Pa;
    StateGrid(row, col) = feasibilityStateCode(ptsBase(i), activeScenario.params);
end

figure;
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
if ~isempty(fieldnames(activeCfg)) && ~isempty(fieldnames(rec))
    title(tl, sprintf(['Constraint and Performance Maps at %.1f mm Pitch\n' ...
        'Comparison-based reference: %s at %.0f V and %.0f L/min'], ...
        baselinePitch, activeCfg.label, rec.voltage_V, rec.totalFlow_Lpm));
elseif ~isempty(fieldnames(rec))
    title(tl, sprintf(['Constraint and Performance Maps at %.1f mm Pitch\n' ...
        'Comparison-based reference at %.0f V and %.0f L/min'], ...
        baselinePitch, rec.voltage_V, rec.totalFlow_Lpm));
else
    title(tl, sprintf('Constraint and Performance Maps at %.1f mm Pitch', ...
        baselinePitch));
end

ax = nexttile;
imagesc(ax, voltVals, flowVals, Qgrid);
set(ax, 'YDir', 'normal');
colorbar(ax);
xlabel(ax, 'Voltage per branch (V)');
ylabel(ax, 'Total flow (L/min)');
title(ax, 'Heat to Bed (W)');
hold(ax, 'on');

ax = nexttile;
imagesc(ax, voltVals, flowVals, CurrentGrid);
set(ax, 'YDir', 'normal');
colorbar(ax);
hold(ax, 'on');
contour(ax, voltVals, flowVals, CurrentGrid, ...
    [activeScenario.params.limits.maxTotalCurrent_A ...
    activeScenario.params.limits.maxTotalCurrent_A], ...
    'k--', 'LineWidth', 1.2);
xlabel(ax, 'Voltage per branch (V)');
ylabel(ax, 'Total flow (L/min)');
title(ax, 'Total Current (A)');

ax = nexttile;
imagesc(ax, voltVals, flowVals, TairGrid);
set(ax, 'YDir', 'normal');
colorbar(ax);
hold(ax, 'on');
contour(ax, voltVals, flowVals, TairGrid, ...
    [activeScenario.params.limits.maxAirOutletTemp_C ...
    activeScenario.params.limits.maxAirOutletTemp_C], ...
    'k--', 'LineWidth', 1.2);
xlabel(ax, 'Voltage per branch (V)');
ylabel(ax, 'Total flow (L/min)');
title(ax, 'Outlet Air Temperature (C)');

ax = nexttile;
imagesc(ax, voltVals, flowVals, WireGrid);
set(ax, 'YDir', 'normal');
colorbar(ax);
hold(ax, 'on');
contour(ax, voltVals, flowVals, WireGrid, ...
    [activeScenario.params.limits.maxWireTemp_C ...
    activeScenario.params.limits.maxWireTemp_C], ...
    'k--', 'LineWidth', 1.2);
xlabel(ax, 'Voltage per branch (V)');
ylabel(ax, 'Total flow (L/min)');
title(ax, 'Max Wire Temperature (C)');

ax = nexttile;
imagesc(ax, voltVals, flowVals, DPgrid);
set(ax, 'YDir', 'normal');
colorbar(ax);
hold(ax, 'on');
contour(ax, voltVals, flowVals, DPgrid, ...
    [activeScenario.params.limits.maxPressureDrop_Pa ...
    activeScenario.params.limits.maxPressureDrop_Pa], ...
    'k--', 'LineWidth', 1.2);
xlabel(ax, 'Voltage per branch (V)');
ylabel(ax, 'Total flow (L/min)');
title(ax, 'Total Pressure Drop (Pa)');

ax = nexttile;
imagesc(ax, voltVals, flowVals, StateGrid);
set(ax, 'YDir', 'normal');
hold(ax, 'on');
stateMap = [0.10 0.60 0.20; 0.80 0.20 0.20; 0.93 0.69 0.13; ...
    0.49 0.18 0.56; 0.00 0.45 0.74; 0.64 0.08 0.18; 0.85 0.33 0.10; 0.50 0.50 0.50];
colormap(ax, stateMap);
clim(ax, [1 8]);
cb = colorbar(ax);
cb.Ticks = 1:8;
cb.TickLabels = {'feasible', 'current', 'air out', 'wire temp', 'pressure', 'water loss', 'hole vel', 'heat shortfall'};
xlabel(ax, 'Voltage per branch (V)');
ylabel(ax, 'Total flow (L/min)');
title(ax, 'Constraint State Map');

if ~isempty(fieldnames(rec))
    recColor = 'w';
    for ax = findall(gcf, 'Type', 'axes').'
        if isa(ax, 'matlab.graphics.axis.Axes') && ~strcmp(ax.Tag, 'Colorbar')
            xline(ax, rec.voltage_V, ':', sprintf('%.0f V', rec.voltage_V), ...
                'Color', recColor, 'LabelVerticalAlignment', 'bottom', ...
                'LabelOrientation', 'horizontal', 'LineWidth', 1.0);
            yline(ax, rec.totalFlow_Lpm, ':', sprintf('%.0f L/min', rec.totalFlow_Lpm), ...
                'Color', recColor, 'LabelHorizontalAlignment', 'left', ...
                'LabelOrientation', 'horizontal', 'LineWidth', 1.0);
            plot(ax, rec.voltage_V, rec.totalFlow_Lpm, 'p', ...
                'MarkerSize', 12, 'MarkerFaceColor', recColor, ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
        end
    end
end

cfg = results.configurationComparisons;
cfgLabels = {cfg.label};
cfgHeat_W = arrayfun(@(s) s.scenario.requiredHeat_W, cfg);
cfgUA_W_K = arrayfun(@(s) s.scenario.lumpedLoss.UA_W_K, cfg);
cfgPower_W = arrayfun(@(s) scenarioReferencePoint(s.scenario).totalPower_W, cfg);
cfgFlow_Lpm = arrayfun(@(s) scenarioReferencePoint(s.scenario).totalFlow_Lpm, cfg);
cfgCurrent_A = arrayfun(@(s) scenarioReferencePoint(s.scenario).totalCurrent_A, cfg);
cfgVoltage_V = arrayfun(@(s) scenarioReferencePoint(s.scenario).voltage_V, cfg);
cfgPitch_mm = arrayfun(@(s) scenarioReferencePoint(s.scenario).pitch_mm, cfg);

figure;
tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
if ~isempty(fieldnames(activeCfg))
    title(tl, sprintf(['Vessel Configuration Comparison\n' ...
        'Downstream outputs use: %s at %.1f mm, %.0f V, %.0f L/min'], ...
        activeCfg.label, rec.pitch_mm, rec.voltage_V, rec.totalFlow_Lpm));
else
    title(tl, 'Vessel Configuration Comparison');
end

ax = nexttile;
b = bar(ax, categorical(cfgLabels), cfgUA_W_K, 'FaceColor', 'flat');
b.CData = vertcat(cfg.color);
ylabel(ax, 'UA (W/K)');
title(ax, 'Overall Ambient Conductance');
labelBarSeries(ax, b, '%.2f');

ax = nexttile;
b = bar(ax, categorical(cfgLabels), cfgHeat_W, 'FaceColor', 'flat');
b.CData = vertcat(cfg.color);
ylabel(ax, 'Required heat (W)');
title(ax, 'Heat Required at Setpoint');
labelBarSeries(ax, b, '%.1f W');

ax = nexttile;
b = bar(ax, categorical(cfgLabels), cfgPower_W, 'FaceColor', 'flat');
b.CData = vertcat(cfg.color);
ylabel(ax, 'Recommended power (W)');
title(ax, 'Comparison-Based Recommended Power');
labelBarSeries(ax, b, '%.1f W');

ax = nexttile;
xPos = 1:numel(cfg);
yyaxis(ax, 'left');
bar(ax, xPos, cfgFlow_Lpm, 'FaceColor', [0.80 0.80 0.80]);
ylabel(ax, 'Recommended flow (L/min)');
for i = 1:numel(cfg)
    text(ax, xPos(i), cfgFlow_Lpm(i), sprintf('%.0f L/min\n%.0f V, %.1f mm', ...
        cfgFlow_Lpm(i), cfgVoltage_V(i), cfgPitch_mm(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 8);
end
yyaxis(ax, 'right');
plot(ax, xPos, cfgCurrent_A, 'ko-', 'LineWidth', 1.4, 'MarkerFaceColor', 'k');
yline(ax, activeScenario.params.limits.maxTotalCurrent_A, 'r--', '15 A cap');
ylabel(ax, 'Recommended current (A)');
for i = 1:numel(cfg)
    text(ax, xPos(i), cfgCurrent_A(i), sprintf('%.2f A', cfgCurrent_A(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 8, 'Color', 'k');
end
set(ax, 'XTick', xPos, 'XTickLabel', cfgLabels);
xlim(ax, [0.5, numel(cfg) + 0.5]);
title(ax, 'Flow, Current, Voltage, and Pitch by Vessel Configuration');

wireStudy = results.wireDiameterSweep;
figure;
tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, ['Wire Diameter Trade Study: Criteria-Based Feasible Selection' newline ...
    'Each curve shows the lexicographically selected feasible point at each wire diameter']);

ax1 = nexttile;
hold(ax1, 'on');
for i = 1:numel(wireStudy)
    plot(ax1, wireStudy(i).diameter_mm, wireStudy(i).optimalPower_W, '-o', ...
        'Color', wireStudy(i).color, 'LineWidth', 1.5, ...
        'MarkerFaceColor', wireStudy(i).color, ...
        'DisplayName', wireStudy(i).label);
    if ~isnan(wireStudy(i).bestIndex)
        idx = wireStudy(i).bestIndex;
        plot(ax1, wireStudy(i).diameter_mm(idx), wireStudy(i).optimalPower_W(idx), ...
            'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'y', 'HandleVisibility', 'off');
    end
end
xlabel(ax1, 'Wire diameter (mm)');
ylabel(ax1, 'Selected total power (W)');
title(ax1, 'Power at Selected Feasible Point');
legend(ax1, 'Location', 'best');
grid(ax1, 'on');

ax2 = nexttile;
hold(ax2, 'on');
for i = 1:numel(wireStudy)
    plot(ax2, wireStudy(i).diameter_mm, wireStudy(i).optimalFlow_Lpm, '-o', ...
        'Color', wireStudy(i).color, 'LineWidth', 1.5, ...
        'MarkerFaceColor', wireStudy(i).color, ...
        'DisplayName', wireStudy(i).label);
    if ~isnan(wireStudy(i).bestIndex)
        idx = wireStudy(i).bestIndex;
        plot(ax2, wireStudy(i).diameter_mm(idx), wireStudy(i).optimalFlow_Lpm(idx), ...
            'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'y', 'HandleVisibility', 'off');
    end
end
xlabel(ax2, 'Wire diameter (mm)');
ylabel(ax2, 'Selected flow rate (L/min)');
title(ax2, 'Flow at Selected Feasible Point');
grid(ax2, 'on');

ax3 = nexttile;
hold(ax3, 'on');
for i = 1:numel(wireStudy)
    plot(ax3, wireStudy(i).diameter_mm, wireStudy(i).optimalCurrent_A, '-o', ...
        'Color', wireStudy(i).color, 'LineWidth', 1.5, ...
        'MarkerFaceColor', wireStudy(i).color, ...
        'DisplayName', wireStudy(i).label);
    if ~isnan(wireStudy(i).bestIndex)
        idx = wireStudy(i).bestIndex;
        plot(ax3, wireStudy(i).diameter_mm(idx), wireStudy(i).optimalCurrent_A(idx), ...
            'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'y', 'HandleVisibility', 'off');
    end
end
xlabel(ax3, 'Wire diameter (mm)');
ylabel(ax3, 'Selected current (A)');
title(ax3, 'Current at Selected Feasible Point');
grid(ax3, 'on');

if ~isempty(fieldnames(activeScenario.recommendedOperatingPoint))
    op = activeScenario.recommendedOperatingPoint;
    heaterX_air_m = linspace(0, activeScenario.params.heaterTube.length_m, ...
        numel(op.branch.heater.airProfile_C));
    heaterX_wire_m = linspace(activeScenario.params.heaterTube.length_m / ...
        numel(op.branch.heater.wireProfile_C), ...
        activeScenario.params.heaterTube.length_m, ...
        numel(op.branch.heater.wireProfile_C));
    aerationX_m = linspace(0, activeScenario.params.aeration.length_m, ...
        numel(op.branch.aeration.airProfile_C));

    figure;
    tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'Recommended-Point Thermal Profiles');

    ax = nexttile;
    plot(ax, heaterX_air_m, op.branch.heater.airProfile_C, 'b-', 'LineWidth', 1.5);
    hold(ax, 'on');
    plot(ax, heaterX_wire_m, op.branch.heater.wireProfile_C, 'r-', 'LineWidth', 1.5);
    xlabel(ax, 'Axial position in heater tube (m)');
    ylabel(ax, 'Temperature (C)');
    title(ax, 'Heater Tube Air and Wire Profiles');
    legend(ax, 'Air', 'Wire', 'Location', 'best');
    grid(ax, 'on');

    ax = nexttile;
    plot(ax, aerationX_m, op.branch.aeration.airProfile_C, 'k-', 'LineWidth', 1.5);
    xlabel(ax, 'Axial position in aeration tube (m)');
    ylabel(ax, 'Air temperature (C)');
    title(ax, 'Aeration Tube Air Cooling Profile');
    grid(ax, 'on');
end

if activeScenario.lumpedLoss.lumpedStrictlyValid || ...
        ~isempty(fieldnames(activeScenario.transient))
    figure;
    nTiles = activeScenario.lumpedLoss.lumpedStrictlyValid + ...
        ~isempty(fieldnames(activeScenario.transient));
    tiledlayout(nTiles, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    if activeScenario.lumpedLoss.lumpedStrictlyValid
        ax = nexttile;
        plot(ax, activeScenario.lumpedLoss.cooldown_time_h, ...
            activeScenario.lumpedLoss.cooldown_C, ...
            'k', 'LineWidth', 1.5);
        hold(ax, 'on');
        yline(ax, activeScenario.params.bin.setpoint_C, 'b--');
        yline(ax, activeScenario.params.environment.greenhouseAir_C, 'r--');
        xlabel(ax, 'Time (h)');
        ylabel(ax, 'Bin temperature (C)');
        title(ax, 'Lumped Cooldown to Greenhouse Air');
        legend(ax, 'Lumped cooldown', 'Setpoint', 'Greenhouse air', ...
            'Location', 'best');
    end

    if ~isempty(fieldnames(activeScenario.transient))
        ax = nexttile;
        plot(ax, activeScenario.transient.time_h, ...
            activeScenario.transient.bottom_C, ...
            'b', 'LineWidth', 1.5);
        hold(ax, 'on');
        plot(ax, activeScenario.transient.time_h, ...
            activeScenario.transient.top_C, ...
            'r', 'LineWidth', 1.5);
        yyaxis(ax, 'right');
        plot(ax, activeScenario.transient.time_h, ...
            activeScenario.transient.voltage_V, ...
            'k--', 'LineWidth', 1.0);
        ylabel(ax, 'Voltage (V)');
        yyaxis(ax, 'left');
        xlabel(ax, 'Time (h)');
        ylabel(ax, 'Bed temperature (C)');
        title(ax, 'Two-Node PID Transient');
        legend(ax, 'Bottom', 'Top', 'Voltage', 'Location', 'best');
    end
end
end

function labelBarSeries(ax, barHandle, fmt)
xVals = barHandle.XEndPoints;
yVals = barHandle.YEndPoints;
yLim = ylim(ax);
yOffset = 0.02 * max(yLim(2) - yLim(1), 1);

for i = 1:numel(yVals)
    text(ax, xVals(i), yVals(i) + yOffset, sprintf(fmt, yVals(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 8);
end
end

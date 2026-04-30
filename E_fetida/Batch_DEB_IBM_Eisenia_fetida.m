%% BATCH POPULATION IBM (DEB-based) ??? Eisenia fetida

% run_Eisenia_fetida.m, mydata_Eisenia_fetida.m,
% pars_init_Eisenia_fetida.m, and predict_Eisenia_fetida.m
% are found here: https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/entries_web/Eisenia_fetida/Eisenia_fetida_res.html
% DEBtool is from https://github.com/add-my-pet/DEBtool_M 
% Adds TOC% and TN% vs time (dry basis) for solids in bin:
%   solids = remaining substrate + accumulated faeces (castings) in current cycle
%
% 
% - Substrate S is wet mass (g).
% - Moisture fixed at 60% => DM = 0.40 (tunable).
% - TOC% and TN% are defined on dry mass basis.
% - Faeces is the solid unassimilated fraction (modeled via kap_P in pars_init_Eisenia_fetida).
% - Assimilated C that is not in solids is assumed to leave as CO2 (mineralization).
% - Assimilated N not in solids is assumed to leave as mineral N (not counted in solids),
%   unless you explicitly model immobilization (not done here).
%
% Operational assumption:
% - Bin contents replaced each cycle, worms left alone.
% - Eggs are removed at reset, incubated externally, and reintroduced upon hatching.

if ~exist('skip_clear', 'var') || ~skip_clear
    clear; clc; close all;
end

%% check DEBtool availability
required_functions = {'parscomp_st', 'tempcorr', 'get_tp', 'get_tm_s', 'reprod_rate'};
missing = {};
for i = 1:length(required_functions)
    if ~exist(required_functions{i}, 'file')
        missing{end+1} = required_functions{i}; %#ok<AGROW>
    end
end
if ~isempty(missing)
    error('Missing required DEBtool functions: %s\nEnsure DEBtool is in MATLAB path.', strjoin(missing, ', '));
end

%% settings
if ~exist('T_C', 'var'); T_C = 20; % C
end      
if ~exist('t_end', 'var'); t_end = 90; % days
end     
dt          = 1;       % days
if ~exist('cycle_len', 'var'); cycle_len = 45; % days
end      

% batch substrate - wet mass (total available food mass per cycle) 
S0          = 121512;    % g wet substrate per cycle [S0/(1215.12g/gal)] = 100 gallons) assuming a bulk density of 0.321kg/L for precomposted (2 months) tomato crop residue and almond shells (75:25, v:v) per http://dx.doi.org/10.1016/j.biortech.2012.05.028
K_S         = 0.15*S0; % g, half-saturation for f(S)=S/(S+K_S)

% Initial population (adults)
N_adult0    = 10000; %60756;     % density of 0.5 worm / g (very high; HEVSTOW had 40% conversion of substrate within 20 days for 50 worms/200g = 0.25 worm/g and 100% conversion at a higher, unspecified worm density; see https://doi.org/10.1007/s13399-021-01522-w)
Ww_adult0   = 0.6;     % g wet weight per initial adult (approx) from run_Eisenia_fetida.m

% optional initial juveniles
N_juv0      = 0;
Ww_juv0     = 0.15;    % g wet weight per initial juvenile

% mortality toggle
use_mortality = true;
use_deb_aging = true;

% egg-only control every cycle to keep worm population stable next cycle
if ~exist('enable_egg_control', 'var'); enable_egg_control = true; end
enable_egg_control_debug = true;
control_horizon_cycles = 3;

% empirical ingestion
use_empirical_ingestion = true;
ingestion_rate_percent = 0.30;  % 25-35% of body weight per day (scaled by f) from https://composting.ces.ncsu.edu/vermicomposting-2/wormy-facts-and-interesting-tidbits/
% suppress_summary: silence console output for fitting
if ~exist('suppress_summary', 'var'); suppress_summary = false; end
% suppress_plots: skip all figure/plot generation during fitting calls
if ~exist('suppress_plots', 'var'); suppress_plots = false; end
% fit_mode: deterministic/simplified path for calibration runs
if ~exist('fit_mode', 'var'); fit_mode = logical(suppress_summary) && logical(suppress_plots); end
if fit_mode
    use_mortality = false;
    use_deb_aging = false;
    enable_egg_control = false;
end
if ~suppress_summary
    fprintf('Using empirical ingestion: %.0f%% body weight/day\n', ingestion_rate_percent*100);
    if fit_mode
        fprintf('Fit mode: yes\n');
    else
        fprintf('Fit mode: no\n');
    end
end

% reproduction correction factor
reprod_correction = 0.36 / 0.888;
if ~suppress_summary
    fprintf('Reproduction correction factor: %.4f\n', reprod_correction);
end

% egg hatch survival fraction within cocoons
% https://www.sciencedirect.com/science/article/abs/pii/S0929139396001205
hatch_surv  = 0.64;

% egg incubation time (cocoon to hatch -> 30-75 days)
% https://composting.ces.ncsu.edu/vermicomposting-2/wormy-facts-and-interesting-tidbits/
egg_incubation_days = 52;

% progress printing
if ~exist('progressEveryDays', 'var'); progressEveryDays = 1; end
rng(1);

%% TOC/TN settings
% given by you: wet mass with moisture 60%
if ~exist('moisture', 'var'); moisture = 0.60; end
if ~exist('moisture_t', 'var'); moisture_t = []; end      % optional time series (days)
if ~exist('moisture_series', 'var'); moisture_series = []; end % optional time series (wet basis, 0-1)
moisture_init = get_moisture(0, moisture, moisture_t, moisture_series);
DM = 1 - moisture_init;        % dry matter fraction = 0.40 (initial)

% set to match Table 3 of http://dx.doi.org/10.1016/j.biortech.2012.05.028
if ~exist('TOC0_dw', 'var'); TOC0_dw = 0.349; end          % 34.9% of dry mass is carbon
if ~exist('TN0_dw', 'var'); TN0_dw  = 0.0270; end         % 2.70% of dry mass is nitrogen

% bin emptied each cycle since substrate replaced entirely
empty_bin_each_cycle = true;

% chemistry knobs (constants; fit once for given conditions)
% k_N_loss: first-order nitrogen loss from solids (NH3/leachate)
% k_C_decay: first-order carbon decay from solids (microbial mineralization)
% k_N_immob: first-order immobilization of mineral N into solids (net gain in solid N)
% set to match Tables 2 & 3 of http://dx.doi.org/10.1016/j.biortech.2009.07.030
k_N_loss  = 1.0e-5; % 1/d (fitted for 0-20 d dataset from)
k_C_decay = 0.00533577;  % 1/d (fitted for 0-20 d dataset)
k_N_immob = 0.0595883;  % 1/d (fitted for 0-20 d dataset)
% Source (DEB flux partitioning and mineral flux context):
% https://doi.org/10.1017/CBO9780511805400
% https://academic.oup.com/conphys/article/4/1/cow023/2951363
% Source (vermicomposting N transformation/mineralization behavior):
% https://pubmed.ncbi.nlm.nih.gov/31096104/
% https://pubmed.ncbi.nlm.nih.gov/23831778/
frac_N_loss_assim = 0.000333021; % fitted for 0-20 d dataset
frac_N_immob_assim = 0.691235;   % fitted for 0-20 d dataset
CN_immob_target = 10;

% optional: parameter modulation by temperature and moisture
% param_coeffs is at least 5x3 [b0 bT bM] for:
% [k_C_decay; k_N_loss; k_N_immob; frac_N_loss_assim; frac_N_immob_assim]
% Optional 4th column bTM adds temperature-moisture interaction:
% sigmoid(b0 + bT*xT + bM*xM + bTM*xT*xM), where xT=(T-refT)/param_T_scale, xM=(M-refM)/param_M_scale.
% Optional rows 6:12 control chemistry extensions:
% [kC_sub_mult; kC_faec_mult; kNloss_sub_mult; kNloss_faec_mult; kNimmob_sub_mult; kNimmob_faec_mult; k_DM_other]
if ~exist('param_coeffs', 'var'); param_coeffs = []; end
if ~exist('param_T_ref', 'var'); param_T_ref = 20; end
if ~exist('param_M_ref', 'var'); param_M_ref = 0.60; end
if ~exist('param_lb', 'var'); param_lb = [1e-4, 1e-5, 1e-5, 0.00, 0.00]; end
if ~exist('param_ub', 'var'); param_ub = [0.20, 0.10, 0.10, 0.50, 0.90]; end
if ~exist('param_T_scale', 'var'); param_T_scale = 1.0; end
if ~exist('param_M_scale', 'var'); param_M_scale = 1.0; end
if ~exist('param_extra_lb', 'var'); param_extra_lb = [0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.00]; end
if ~exist('param_extra_ub', 'var'); param_extra_ub = [2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 0.08]; end

% chemistry extensions (defaults if not estimated via param_coeffs extra rows)
chem_extra_defaults = struct();
chem_extra_defaults.kC_sub_mult = 1.0;
chem_extra_defaults.kC_faec_mult = 1.0;
chem_extra_defaults.kNloss_sub_mult = 1.0;
chem_extra_defaults.kNloss_faec_mult = 1.0;
chem_extra_defaults.kNimmob_sub_mult = 1.0;
chem_extra_defaults.kNimmob_faec_mult = 1.0;
chem_extra_defaults.k_DM_other = 0.0;
% Source (density dependence in vermicomposting outcomes; used to motivate rho terms):
% https://pubmed.ncbi.nlm.nih.gov/27023819/
% https://www.nature.com/articles/s41598-022-17855-z
use_density_feedback = true;
rho_half = 0.5;
k_rho_f = 1.5;
kC_rho = 1.0;
kC_act = 1.0;
kNloss_rho = 0.8;
kNloss_act = 1.2;
kNimmob_rho = 0.3;
kNimmob_act = 0.2;
sel_C0 = 0.05;
sel_N0 = 0.20;
sel_C_rho = 0.10;
sel_N_rho = 0.20;
sel_C_act = 0.10;
sel_N_act = 0.25;
if fit_mode
    use_density_feedback = false;
end

%% food limitation settings
% thresholds for "food limitation" reporting:
% - catabolism trigger (growth stops or reverses): f <= l_T
% - ingestion too low: per-worm ingestion < ingestion_floor_frac of max-per-worm rate
ingestion_floor_frac = 0.20;

%% load DEB parameters
[data, auxData, metaData] = mydata_Eisenia_fetida; %#ok<ASGLU>
[par, metaPar] = pars_init_Eisenia_fetida(metaData); %#ok<ASGLU>
cPar = parscomp_st(par);

% temperature correction
TC = tempcorr(C2K(T_C), par.T_ref, par.T_A);

% pull parameters (from pars_init_Eisenia_fetida.m)
k_M   = cPar.k_M;     % 1/d
L_m   = cPar.L_m;     % cm
l_T   = cPar.l_T;     % -
v     = par.v;        % cm/d
p_Am  = cPar.p_Am;    % J/d/cm^2

% temperature-corrected rate parameters (Arrhenius)
k_M_T = k_M * TC;
v_T   = v * TC;
k_J_T = par.k_J * TC;
p_Am_T = p_Am * TC;

kap_X = par.kap_X;    % digestion efficiency food -> reserve
kap_P = par.kap_P;    % faecation efficiency food -> faeces

% mu_X and w_X (used in mechanistic ingestion branch)
if isfield(par,'mu_X'); mu_X = par.mu_X;
elseif isfield(cPar,'mu_X'); mu_X = cPar.mu_X;
else; error('mu_X not found in par or cPar.'); end

if isfield(par,'w_X'); w_X = par.w_X;
elseif isfield(cPar,'w_X'); w_X = cPar.w_X;
else; error('w_X not found in par or cPar.'); end

w = cPar.w;

% life history at reference f=1 (used for birth size + puberty size)
f_ref = 1;

pars_tp = [cPar.g cPar.k l_T cPar.v_Hb cPar.v_Hp];
[t_p, t_b, l_p, l_b, info_tp] = get_tp(pars_tp, f_ref);
if info_tp ~= 1
    error('get_tp failed (info=%g).', info_tp);
end
L_b = L_m * l_b;
L_p = L_m * l_p;

% egg incubation time override
a_b_T = egg_incubation_days;

% lifespan estimate -> crude constant hazard
pars_tm = [cPar.g; l_T; par.h_a / k_M^2; par.s_G];
t_m = get_tm_s(pars_tm, f_ref, l_b);
a_m_T = max(t_m / k_M_T, 10);
mu_adult = 1 / a_m_T;
mu_juv   = 1 / max(0.7*a_m_T, 7);
mu_egg   = 1 / max(0.2*a_m_T, 2);
h_a_T = par.h_a * TC^2;
s_G = par.s_G;

% reproduction parameters
pars_R = [par.kap; par.kap_R; cPar.g; k_J_T; k_M_T; cPar.L_T; v_T; cPar.U_Hb; cPar.U_Hp];

% check whether reprod_rate supports vector inputs
use_vector_reprod = false;
if ~fit_mode
    try
        test_L = [L_b; L_p];
        test_R = reprod_rate(test_L, f_ref, pars_R);
        if isvector(test_R) && numel(test_R) == numel(test_L)
            use_vector_reprod = true;
        end
    catch
        use_vector_reprod = false;
    end
end

%% initialize population
L_adult = (Ww_adult0 ./ (1 + f_ref*w)).^(1/3) * ones(N_adult0,1);
L_juv   = (Ww_juv0   ./ (1 + f_ref*w)).^(1/3) * ones(N_juv0,1);
q_adult = zeros(N_adult0,1);
h_adult = zeros(N_adult0,1);
q_juv   = zeros(N_juv0,1);
h_juv   = zeros(N_juv0,1);

% eggs as age bins:
% - bin eggs (in substrate during cycle)
% - incubator eggs (external incubation across cycles)
nEggBins = ceil((max(cycle_len, a_b_T) + 5) / dt) + 5;
egg_bins_bin = zeros(nEggBins,1);
egg_bins_inc = zeros(nEggBins,1);

% batch substrate remaining (wet mass)
S = S0;
next_reset = cycle_len;

%% TOC/TN state per cycle
% substrate solids pool (dry mass + C/N)
DM_sub = (S0 * DM);
C_sub  = DM_sub * TOC0_dw;
N_sub  = DM_sub * TN0_dw;

% faeces solids pool (dry mass + C/N), plus wet mass for plotting
DM_faec = 0;
C_faec  = 0;
N_faec  = 0;
P_wet   = 0;  % g wet

%% output arrays
t = (0:dt:t_end).';
nt = numel(t);

J_ing   = zeros(nt,1);  % g/d wet ingestion
J_faec  = zeros(nt,1);  % g/d wet faeces
F_cum   = zeros(nt,1);  % g cumulative wet faeces
G_rate  = zeros(nt,1);  % g/d total wet biomass growth
Negg    = zeros(nt,1);
Negg_bin = zeros(nt,1);
Negg_inc = zeros(nt,1);
Njuv    = zeros(nt,1);
Nad     = zeros(nt,1);
S_rem   = zeros(nt,1);
f_hist  = zeros(nt,1);
keep_frac = NaN(nt,1);
keep_eggs = NaN(nt,1);
keep_worms = NaN(nt,1);
expected_worms_next = NaN(nt,1);
expected_worms_survivors = NaN(nt,1);
expected_hatch_all = NaN(nt,1);
last_expected_worms_next = NaN;
last_expected_worms_survivors = NaN;
last_expected_hatch_all = NaN;
last_keep_frac = NaN;
last_keep_eggs = NaN;
last_keep_worms = NaN;
food_limit_flag = false;
food_limit_time = NaN;
food_limit_worms = NaN;
food_limit_density = NaN;

% TOC/TN outputs (dry basis, for solids present in bin in current cycle)
TOC_pct = zeros(nt,1);
TN_pct  = zeros(nt,1);
TOC_pct_wet = zeros(nt,1);
TN_pct_wet  = zeros(nt,1);
CNRatio = zeros(nt,1);

Wtot_prev = sum(L_adult.^3) * (1 + f_ref*w) + sum(L_juv.^3) * (1 + f_ref*w);

progressEveryK = max(1, round(progressEveryDays/dt));
tic;

%% main
for k = 1:nt
    tk = t(k);
    moisture_now = get_moisture(tk, moisture, moisture_t, moisture_series);
    DM_now = 1 - moisture_now;

    % batch event: egg control + substrate regenerated (bin emptied)
    if tk >= next_reset - 1e-12
        % egg-only control to keep worm population stable next cycle
        [egg_bins_bin, egg_bins_inc, L_adult, L_juv, keep_fraction, keep_eggs(k), keep_worms(k), ...
            expected_worms_survivors(k), expected_hatch_all(k), expected_worms_next(k)] = ...
            apply_cycle_egg_control(egg_bins_bin, egg_bins_inc, L_adult, L_juv, enable_egg_control, ...
            cycle_len, dt, a_b_T, mu_egg, hatch_surv, mu_adult, mu_juv, control_horizon_cycles, enable_egg_control_debug, tk, fit_mode);
        keep_frac(k) = keep_fraction;
        last_keep_frac = keep_fraction;
        last_keep_eggs = keep_eggs(k);
        last_keep_worms = keep_worms(k);
        last_expected_worms_survivors = expected_worms_survivors(k);
        last_expected_hatch_all = expected_hatch_all(k);
        last_expected_worms_next = expected_worms_next(k);

        % replace substrate entirely
        S = S0;

        % reset solids chemistry to fresh substrate
        DM_sub = (S0 * DM);
        C_sub  = DM_sub * TOC0_dw;
        N_sub  = DM_sub * TN0_dw;
        DM_faec = 0;
        C_faec  = 0;
        N_faec  = 0;

        if empty_bin_each_cycle
            P_wet  = 0;
        end

        next_reset = next_reset + cycle_len;
    end

    % functional response
    % Source: DEB functional response f = X/(X+X_K), with density modifier added here
    % to reflect crowding/interference in vermicomposting systems.
    % https://academic.oup.com/conphys/article/4/1/cow023/2951363
    % https://pubmed.ncbi.nlm.nih.gov/27023819/
    total_worms_now = numel(L_adult) + numel(L_juv);
    rho_w = total_worms_now / max(S + P_wet, 1);
    rho_idx = rho_w / (rho_w + rho_half);
    K_S_eff = K_S;
    if use_density_feedback
        K_S_eff = K_S * (1 + k_rho_f * rho_idx);
    end
    f = S / (S + K_S_eff);
    f = min(max(f,0),1);

    % eggs: age forward (shift bins), mortality, hatch (incubator only)
    % NUMERICAL METHOD: discrete age-bin shift (deterministic aging)
    egg_bins_bin(2:end) = egg_bins_bin(1:end-1);
    egg_bins_bin(1) = 0;
    egg_bins_inc(2:end) = egg_bins_inc(1:end-1);
    egg_bins_inc(1) = 0;

    if use_mortality
        if fit_mode
            % Deterministic expectation in fitting mode.
            egg_surv = exp(-mu_egg*dt);
            egg_bins_bin = egg_bins_bin * egg_surv;
            egg_bins_inc = egg_bins_inc * egg_surv;
        else
            % NUMERICAL METHOD: stochastic survival via binomial thinning
            egg_bins_bin = binornd(round(egg_bins_bin), exp(-mu_egg*dt));
            egg_bins_inc = binornd(round(egg_bins_inc), exp(-mu_egg*dt));
        end
    end

    hatch_bin = floor(a_b_T/dt) + 1;
    if hatch_bin <= nEggBins
        n_hatch = egg_bins_inc(hatch_bin);
        egg_bins_inc(hatch_bin) = 0;
        if n_hatch > 0
            if fit_mode
                n_new = round(max(n_hatch, 0) * hatch_surv);
            else
                n_new = binornd(round(n_hatch), hatch_surv);
            end
            if n_new > 0
                L_juv = [L_juv; L_b * ones(n_new,1)];
                q_juv = [q_juv; zeros(n_new,1)];
                h_juv = [h_juv; zeros(n_new,1)];
            end
        end
    end

    % juveniles: grow + ingest + faeces + mortality + mature
    JX_juv = 0; JP_juv = 0;

    % precompute growth constants for this step
    L_i_common = L_m * max(f - l_T, 0);
    ir_B = 3/k_M_T + 3 * f * L_m / v_T;
    r_B  = 1 / ir_B;

    if ~isempty(L_juv)
        if ~isempty(L_juv)
            L_i = max(L_i_common, L_p * 1.01);
            dLdt = r_B .* (L_i - L_juv);
            % NUMERICAL METHOD: forward Euler update for growth ODE
            if use_mortality && use_deb_aging
                [q_juv, h_juv] = deb_aging_step(L_juv, q_juv, h_juv, dLdt, f, L_m, v_T, h_a_T, s_G, dt);
                % NUMERICAL METHOD: stochastic survival using boolean mask (Bernoulli trials per individual - we kill worms off here)
                
                % Eliminates each worm by eliminating the assigned
                % state vector of each descriptive variable

                % survival probability is exp(-integral(h_juv*dt)), 
                % but given the miniscule timestep, it is simplified to
                % exp(-h_juv*dt) for simplicity (see S(t) in
                % https://web.stanford.edu/class/stats305b/notes/Survival_I.html)
                
                p_survive = exp(-h_juv*dt);
                if fit_mode
                    survive = deterministic_survival_mask(p_survive);
                else
                    survive = rand(size(L_juv)) < p_survive;
                end
                L_juv = L_juv(survive);
                q_juv = q_juv(survive);
                h_juv = h_juv(survive);
                dLdt = dLdt(survive);
            elseif use_mortality
                % NUMERICAL METHOD: stochastic survival (Bernoulli trials per individual)
                p_survive = exp(-mu_juv*dt) * ones(size(L_juv));
                if fit_mode
                    survive = deterministic_survival_mask(p_survive);
                else
                    survive = rand(size(L_juv)) < p_survive;
                end
                L_juv = L_juv(survive);
                dLdt = dLdt(survive);
            end

            L_juv = max(L_juv + dLdt*dt, 1e-9);

            if use_empirical_ingestion
                Ww_juv = L_juv.^3 .* (1 + f*w);
                J_X = f .* ingestion_rate_percent .* Ww_juv; % g wet/d per individual
                JX_juv = sum(J_X);
            else
                p_A = f .* p_Am_T .* (L_juv.^2);
                J_X_mol = p_A ./ (kap_X .* mu_X);      % mol/d food
                J_X_dry = J_X_mol .* w_X;              % g dry/d food
                J_X = J_X_dry ./ max(DM_now, 1e-12);   % g wet/d food
                JX_juv = sum(J_X);
            end

            J_P = kap_P .* J_X;
            JP_juv = sum(J_P);

            to_adult = L_juv >= L_p;
            if any(to_adult)
                L_adult = [L_adult; L_juv(to_adult)];
                q_adult = [q_adult; q_juv(to_adult)];
                h_adult = [h_adult; h_juv(to_adult)];
                L_juv = L_juv(~to_adult);
                q_juv = q_juv(~to_adult);
                h_juv = h_juv(~to_adult);
            end
        end
    end

    % adults: grow + ingest + faeces + reproduce + mortality
    JX_ad = 0; JP_ad = 0;

    if ~isempty(L_adult)
        if ~isempty(L_adult)
            L_i = L_i_common;
            dLdt = r_B .* (L_i - L_adult);
            % NUMERICAL METHOD: forward Euler update for growth ODE
            if use_mortality && use_deb_aging
                [q_adult, h_adult] = deb_aging_step(L_adult, q_adult, h_adult, dLdt, f, L_m, v_T, h_a_T, s_G, dt);
                % NUMERICAL METHOD: stochastic survival (Bernoulli trials per individual)
                p_survive = exp(-h_adult*dt);
                if fit_mode
                    survive = deterministic_survival_mask(p_survive);
                else
                    survive = rand(size(L_adult)) < p_survive;
                end
                L_adult = L_adult(survive);
                q_adult = q_adult(survive);
                h_adult = h_adult(survive);
                dLdt = dLdt(survive);
            elseif use_mortality
                % NUMERICAL METHOD: stochastic survival (Bernoulli trials per individual)
                p_survive = exp(-mu_adult*dt) * ones(size(L_adult));
                if fit_mode
                    survive = deterministic_survival_mask(p_survive);
                else
                    survive = rand(size(L_adult)) < p_survive;
                end
                L_adult = L_adult(survive);
                dLdt = dLdt(survive);
            end

            L_adult = max(L_adult + dLdt*dt, 1e-9);

            if use_empirical_ingestion
                Ww_adult = L_adult.^3 .* (1 + f*w);
                J_X = f .* ingestion_rate_percent .* Ww_adult;
                JX_ad = sum(J_X);
            else
                p_A = f .* p_Am_T .* (L_adult.^2);
                J_X_mol = p_A ./ (kap_X .* mu_X);      % mol/d food
                J_X_dry = J_X_mol .* w_X;              % g dry/d food
                J_X = J_X_dry ./ max(DM_now, 1e-12);   % g wet/d food
                JX_ad = sum(J_X);
            end

            J_P = kap_P .* J_X;
            JP_ad = sum(J_P);

            if fit_mode
                % Reproduction is bypassed in fit mode to avoid DEB puberty
                % feasibility warnings and keep calibration objective smooth.
                eggs_expected = 0;
                n_new_eggs = 0;
            else
                if use_vector_reprod
                    R = reprod_correction .* reprod_rate(L_adult, f, pars_R);
                else
                    R = reprod_correction .* arrayfun(@(Li) reprod_rate(Li, f, pars_R), L_adult);
                end
                eggs_expected = max(sum(R) * dt, 0);
                % NUMERICAL METHOD: stochastic rounding of expected count
                n_new_eggs = floor(eggs_expected) + (rand < (eggs_expected - floor(eggs_expected)));
            end
            if n_new_eggs > 0
                egg_bins_bin(1) = egg_bins_bin(1) + n_new_eggs;
            end
        end
    end

    % substrate update (wet mass)
    JX_tot = JX_juv + JX_ad;      % g wet/d requested
    S_before = S;

    % NUMERICAL METHOD: explicit mass balance with hard cap (algebraic)
    dS_ing_wet = min(JX_tot * dt, S_before);  % cannot ingest more than available
    S = max(S_before - dS_ing_wet, 0);

    % faeces wet production this step
    dF_wet = (JP_juv + JP_ad) * dt;      % g wet faeces produced

    % ingestion/faeces outputs
    J_ing(k)  = dS_ing_wet / dt;         % g/d realized
    J_faec(k) = dF_wet / dt;             % g/d

    if k == 1
        F_cum(k) = dF_wet;
    else
        F_cum(k) = F_cum(k-1) + dF_wet;
    end

    % growth output
    Wtot = sum(L_adult.^3) * (1 + f*w) + sum(L_juv.^3) * (1 + f*w);
    G_rate(k) = (Wtot - Wtot_prev)/dt;
    Wtot_prev = Wtot;

    % food-limitation detection (first occurrence)
    if ~food_limit_flag
        total_worms_now = numel(L_adult) + numel(L_juv);
        per_worm_ing = J_ing(k) / max(total_worms_now, 1);
        Ww_mean = Wtot / max(total_worms_now, 1);
        per_worm_max = ingestion_rate_percent * Ww_mean;
        ingestion_too_low = (per_worm_max > 0) && (per_worm_ing < ingestion_floor_frac * per_worm_max);
        catabolism = (f <= l_T);
        if catabolism || ingestion_too_low
            food_limit_flag = true;
            food_limit_time = tk;
            food_limit_worms = total_worms_now;
            food_limit_density = total_worms_now / max(S + P_wet, 1); % worms per g total wet solids
        end
    end

    % TOC/TN total solids (book-keeping; DEB-consistent?)
    dS_ing_dry = dS_ing_wet * DM_now;

    % remove ingested dry mass from substrate pool
    if DM_sub > 0
        Cfrac_sub = C_sub / DM_sub;
        Nfrac_sub = N_sub / DM_sub;
    else
        Cfrac_sub = TOC0_dw;
        Nfrac_sub = TN0_dw;
    end

    C_removed = Cfrac_sub * dS_ing_dry;
    N_removed = Nfrac_sub * dS_ing_dry;

    DM_sub = max(DM_sub - dS_ing_dry, 0);
    C_sub  = max(C_sub - C_removed, 0);
    N_sub  = max(N_sub - N_removed, 0);

    % Selective C/N routing to faeces (empirical extension): earthworms and associated
    % microbes alter C/N composition during gut processing and vermicast formation.
    % https://pubmed.ncbi.nlm.nih.gov/23831778/
    % https://pubmed.ncbi.nlm.nih.gov/31096104/
    act_idx = min(max((JX_tot / max(Wtot, 1e-9)) / max(ingestion_rate_percent, 1e-9), 0), 1);
    if fit_mode
        sel_C = min(max(sel_C0, 0), 0.95);
        sel_N = min(max(sel_N0, 0), 0.95);
    else
        sel_C = min(max(sel_C0 + sel_C_rho * rho_idx + sel_C_act * act_idx, 0), 0.95);
        sel_N = min(max(sel_N0 + sel_N_rho * rho_idx + sel_N_act * act_idx, 0), 0.95);
    end
    C_to_faec = kap_P * (1 - sel_C) * C_removed;
    N_to_faec = kap_P * (1 - sel_N) * N_removed;

    % add faeces solids (unassimilated fraction)
    DM_faec = DM_faec + kap_P * dS_ing_dry;
    C_faec  = C_faec  + C_to_faec;
    N_faec  = N_faec  + N_to_faec;
    P_wet   = P_wet   + dF_wet;

    % first-order decay losses from total solids
    % NUMERICAL METHOD: forward Euler for first-order decay/immobilization
    C_tot = C_sub + C_faec;
    N_tot = N_sub + N_faec;
    DM_tot = max(DM_sub + DM_faec, 1e-9);

    chem_extra = chem_extra_defaults;
    if ~isempty(param_coeffs)
        [k_C_decay, k_N_loss, k_N_immob, frac_N_loss_assim, frac_N_immob_assim, chem_extra] = ...
            eval_param_coeffs(param_coeffs, param_lb, param_ub, T_C, moisture_now, param_T_ref, param_M_ref, ...
            param_T_scale, param_M_scale, param_extra_lb, param_extra_ub, chem_extra_defaults);
    end

    % Activity- and density-dependent turnover rates (empirical extension).
    if fit_mode
        k_C_core = k_C_decay;
        k_N_loss_core = k_N_loss;
        k_N_immob_core = k_N_immob;
    else
        k_C_core = k_C_decay * (1 + kC_rho * rho_idx + kC_act * act_idx);
        k_N_loss_core = k_N_loss * (1 + kNloss_rho * rho_idx + kNloss_act * act_idx);
        k_N_immob_core = k_N_immob * (1 + kNimmob_rho * rho_idx + kNimmob_act * act_idx);
    end

    % Split substrate/faeces chemistry rates (improves identifiability by allowing pool-specific turnover).
    k_C_sub_eff = max(k_C_core * chem_extra.kC_sub_mult, 0);
    k_C_faec_eff = max(k_C_core * chem_extra.kC_faec_mult, 0);
    k_N_loss_sub_eff = max(k_N_loss_core * chem_extra.kNloss_sub_mult, 0);
    k_N_loss_faec_eff = max(k_N_loss_core * chem_extra.kNloss_faec_mult, 0);
    k_N_immob_sub_eff = max(k_N_immob_core * chem_extra.kNimmob_sub_mult, 0);
    k_N_immob_faec_eff = max(k_N_immob_core * chem_extra.kNimmob_faec_mult, 0);

    C_decay_sub = min(k_C_sub_eff * C_sub * dt, C_sub);
    C_decay_faec = min(k_C_faec_eff * C_faec * dt, C_faec);
    C_decay = C_decay_sub + C_decay_faec;

    % DEB-consistent assimilated-N partition (empirical fractions).
    N_assim = max(N_removed - N_to_faec, 0);
    N_loss_sub = min(k_N_loss_sub_eff * N_sub * dt, N_sub);
    N_loss_faec = min(k_N_loss_faec_eff * N_faec * dt, N_faec);
    N_loss_assim = max(frac_N_loss_assim, 0) * N_assim;
    if N_tot > 0
        N_loss_sub = N_loss_sub + N_loss_assim * (N_sub / N_tot);
        N_loss_faec = N_loss_faec + N_loss_assim * (N_faec / N_tot);
    else
        N_loss_sub = N_loss_sub + 0.5 * N_loss_assim;
        N_loss_faec = N_loss_faec + 0.5 * N_loss_assim;
    end
    N_loss_sub = min(max(N_loss_sub, 0), N_sub);
    N_loss_faec = min(max(N_loss_faec, 0), N_faec);

    % Immobilization coupled to C turnover using target microbial C:N.
    N_immob_sub = k_N_immob_sub_eff * (C_decay_sub / max(CN_immob_target, 1e-9));
    N_immob_faec = k_N_immob_faec_eff * (C_decay_faec / max(CN_immob_target, 1e-9));
    N_immob_assim = max(frac_N_immob_assim, 0) * N_assim;
    if C_decay > 0
        N_immob_sub = N_immob_sub + N_immob_assim * (C_decay_sub / C_decay);
        N_immob_faec = N_immob_faec + N_immob_assim * (C_decay_faec / C_decay);
    else
        N_immob_sub = N_immob_sub + 0.5 * N_immob_assim;
        N_immob_faec = N_immob_faec + 0.5 * N_immob_assim;
    end
    N_immob_sub = max(N_immob_sub, 0);
    N_immob_faec = max(N_immob_faec, 0);

    C_sub = max(C_sub - C_decay_sub, 0);
    C_faec = max(C_faec - C_decay_faec, 0);
    N_sub = max(N_sub - N_loss_sub + N_immob_sub, 0);
    N_faec = max(N_faec - N_loss_faec + N_immob_faec, 0);

    N_decay = N_loss_sub + N_loss_faec;
    N_immob = N_immob_sub + N_immob_faec;

    DM_other = max(DM_tot - C_tot - N_tot, 0);
    k_DM_other_eff = max(chem_extra.k_DM_other, 0);
    DM_other_loss = min(k_DM_other_eff * DM_other * dt, DM_other);
    DM_tot_new = max(DM_tot - (C_decay + N_decay - N_immob) - DM_other_loss, 1e-9);
    scaleDM = DM_tot_new / DM_tot;
    DM_sub = max(DM_sub * scaleDM, 0);
    DM_faec = max(DM_faec * scaleDM, 0);

    TOC_pct(k) = 100 * (C_sub + C_faec) / max((DM_sub + DM_faec), 1e-9);
    TN_pct(k)  = 100 * (N_sub + N_faec) / max((DM_sub + DM_faec), 1e-9);
    CNRatio(k) = (C_sub + C_faec) / max((N_sub + N_faec), 1e-12);
    TOC_pct_wet(k) = TOC_pct(k) * (1 - moisture_now);
    TN_pct_wet(k)  = TN_pct(k)  * (1 - moisture_now);

    % state outputs
    Negg_bin(k) = sum(egg_bins_bin);
    Negg_inc(k) = sum(egg_bins_inc);
    Negg(k)  = Negg_bin(k) + Negg_inc(k);
    Njuv(k)  = numel(L_juv);
    Nad(k)   = numel(L_adult);
    S_rem(k) = S;
    f_hist(k)= f;
    if isnan(keep_frac(k))
        keep_frac(k) = last_keep_frac;
        keep_eggs(k) = last_keep_eggs;
        keep_worms(k) = last_keep_worms;
    end
    if isnan(expected_worms_next(k))
        expected_worms_survivors(k) = last_expected_worms_survivors;
        expected_hatch_all(k) = last_expected_hatch_all;
        expected_worms_next(k) = last_expected_worms_next;
    end

    % progress print
    if mod(k, progressEveryK) == 0
        total_worms = Nad(k) + Njuv(k);
        per_worm_ing = J_ing(k) / max(total_worms, 1);
        fprintf('t=%6.1f d | Nad=%6d Njuv=%6d Negg=%8d | S=%8.1f g | JX=%7.2f g/d (%.5f g/d/worm) | TOC=%.2f%% TN=%.2f%% | elapsed=%6.1f s\n', ...
            tk, Nad(k), Njuv(k), Negg(k), S, J_ing(k), per_worm_ing, TOC_pct(k), TN_pct(k), toc);
    end
end

%% plots
if ~suppress_plots
    figure('Color','w','Name','Batch vermicomposter IBM (DEB-based)');

    subplot(2,2,1);
    plot(t, J_faec, 'LineWidth', 1.5); grid on;
    xlabel('time (d)'); ylabel('faeces rate (g/d, wet)');
    title('Faeces production rate vs time');

    subplot(2,2,2);
    plot(t, F_cum, 'LineWidth', 1.5); grid on;
    xlabel('time (d)'); ylabel('cumulative faeces (g, wet)');
    title('Cumulative faeces vs time');

    subplot(2,2,3);
    plot(t, J_ing, 'LineWidth', 1.5); grid on;
    xlabel('time (d)'); ylabel('ingestion rate (g/d, wet)');
    title('Ingestion rate vs time');

    subplot(2,2,4);
    plot(t, G_rate, 'LineWidth', 1.5); grid on;
    xlabel('time (d)'); ylabel('growth rate dW/dt (g/d, wet)');
    title('Growth rate vs time (total wet biomass)');

    figure('Color','w','Name','Population & batch state');
    yyaxis left;
    plot(t, Nad, 'LineWidth', 1.5); hold on;
    plot(t, Njuv,'LineWidth', 1.5);
    plot(t, Negg_bin,'LineWidth', 1.5);
    plot(t, Negg_inc,'LineWidth', 1.5);
    grid on; ylabel('count');
    yyaxis right;
    plot(t, S_rem,'LineWidth', 1.5);
    ylabel('substrate remaining (g wet)');
    xlabel('time (d)');
    title('Population stages & substrate (replaced every 45 days)');
    legend('Adults','Juveniles','Eggs (bin)','Eggs (incubator)','Substrate','Location','best');

    figure('Color','w','Name','Feeding limitation');
    plot(t, f_hist,'LineWidth',1.5); grid on;
    xlabel('time (d)'); ylabel('f');
    title('f(t)=S/(S+K_S), resets with substrate regeneration');

    figure('Color','w','Name','Egg control each cycle');
    yyaxis left;
    plot(t, (Nad + Njuv),'LineWidth',1.5); hold on;
    plot(t, keep_eggs,'LineWidth',1.5);
    ylabel('count');
    yyaxis right;
    plot(t, keep_frac,'LineWidth',1.5);
    ylabel('keep fraction');
    ylim([0 0.5]);
    grid on;
    xlabel('time (d)');
    legend('Worms','Eggs incubating','Keep fraction','Location','best');
    title('Egg control each cycle');

    figure('Color','w','Name','Egg-control prediction (horizon)');
    plot(t, expected_worms_next,'LineWidth',1.5); hold on;
    plot(t, expected_worms_survivors,'LineWidth',1.5);
    plot(t, expected_hatch_all,'LineWidth',1.5);
    grid on;
    xlabel('time (d)');
    ylabel('count');
    legend('Expected worms next cycle','Expected survivors','Expected hatch (controlled)','Location','best');
    title(sprintf('Egg control targets (at cycle boundaries), horizon = %d cycles', control_horizon_cycles));

    figure('Color','w','Name','TOC% and TN% of solids in bin (dry basis)');
    yyaxis left;
    plot(t, TOC_pct,'LineWidth',1.5); grid on;
    ylabel('TOC (% dry mass)');
    yyaxis right;
    plot(t, TN_pct,'LineWidth',1.5);
    ylabel('TN (% dry mass)');
    xlabel('time (d)');
    title('TOC% and TN% (solids = remaining substrate + faeces, current cycle)');
    legend('TOC%','TN%','Location','best');

    figure('Color','w','Name','C/N ratio of solids');
    plot(t, CNRatio,'LineWidth',1.5); grid on;
    xlabel('time (d)'); ylabel('C/N (mass ratio)');
    title('C/N ratio vs time (solids = remaining substrate + faeces)');
end

%% quick summary
% average TOC/TN at the end of each cycle (use last step before reset)
cycle_end_idx = find(abs(mod(t, cycle_len)) < 1e-9 & t > 0) - 1;
cycle_end_idx = cycle_end_idx(cycle_end_idx >= 1);
avg_TOC_end = NaN;
avg_TN_end = NaN;
if ~isempty(cycle_end_idx)
    avg_TOC_end = mean(TOC_pct(cycle_end_idx));
    avg_TN_end = mean(TN_pct(cycle_end_idx));
end

if ~suppress_summary
    fprintf('\n=== Initial state ===\n');
    fprintf('Adults: %d | Juveniles: %d | Eggs: %d\n', Nad(1), Njuv(1), Negg(1));
    fprintf('Cumulative faeces (wet): %.2f g\n', F_cum(1));
    fprintf('Mean ingestion rate (wet): %.3f g/d\n', J_ing(1));
    fprintf('Mean faeces rate (wet): %.3f g/d\n', J_faec(1));
    fprintf('Initial TOC: %.2f%% (dw), Initial TN: %.2f%% (dw)\n\n', TOC_pct(1), TN_pct(1));

    fprintf('=== Final state ===\n');
    fprintf('Adults: %d | Juveniles: %d | Eggs: %d\n', Nad(end), Njuv(end), Negg(end));
    fprintf('Cumulative faeces (wet): %.2f g\n', F_cum(end));
    fprintf('Mean ingestion rate (wet): %.3f g/d\n', mean(J_ing));
    fprintf('Mean faeces rate (wet): %.3f g/d\n', mean(J_faec));
    fprintf('Final TOC: %.2f%% (dw), Final TN: %.2f%% (dw)\n\n', TOC_pct(end), TN_pct(end));
    fprintf('Avg end-of-cycle TOC: %.2f%% (dw), Avg end-of-cycle TN: %.2f%% (dw)\n\n', avg_TOC_end, avg_TN_end);
    if food_limit_flag
        fprintf('Food limitation first detected at t=%.1f d | worms=%d | density=%.4f worms/g total wet solids\n', ...
            food_limit_time, food_limit_worms, food_limit_density);
    else
        fprintf('Food limitation not detected within simulation window.\n');
    end
end

%% local functions
function [egg_bins_bin, egg_bins_inc, L_adult, L_juv, keep_fraction, keep_eggs, keep_worms, ...
    expected_survivors, expected_hatch_all, expected_worms_next] = ...
    apply_cycle_egg_control(egg_bins_bin, egg_bins_inc, L_adult, L_juv, enable_egg_control, ...
    cycle_len, dt, a_b_T, mu_egg, hatch_surv, mu_adult, mu_juv, control_horizon_cycles, enable_debug, tk, deterministic_mode)
% apply_cycle_egg_control: thin eggs only to keep worm population stable next cycle.

total_eggs_bin = sum(egg_bins_bin);
total_eggs_inc = sum(egg_bins_inc);
total_worms = numel(L_adult) + numel(L_juv);
keep_fraction = 1;
expected_survivors = NaN;
expected_hatch_all = NaN;
expected_worms_next = NaN;

if enable_egg_control && (total_eggs_bin > 0 || total_eggs_inc > 0 || total_worms > 0)
    horizon_len = cycle_len * control_horizon_cycles;

    % expected worm survivors to horizon (counts only)
    expected_survivors = numel(L_adult) * exp(-mu_adult * horizon_len) + ...
        numel(L_juv) * exp(-mu_juv * horizon_len);

    % expected hatch within the horizon: incubator (fixed) + bin eggs (controllable)
    nEggBins = numel(egg_bins_bin);
    expected_hatch_inc = 0;
    expected_hatch_bin_all = 0;
    for i = 1:nEggBins
        age = (i-1) * dt;
        time_to_hatch = a_b_T - age;
        if time_to_hatch >= 0 && time_to_hatch <= horizon_len
            expected_hatch_inc = expected_hatch_inc + egg_bins_inc(i) * exp(-mu_egg * time_to_hatch) * hatch_surv;
            expected_hatch_bin_all = expected_hatch_bin_all + egg_bins_bin(i) * exp(-mu_egg * time_to_hatch) * hatch_surv;
        end
    end

    needed_from_new = max(total_worms - expected_survivors - expected_hatch_inc, 0);
    if expected_hatch_bin_all > 0
        keep_fraction = min(1, needed_from_new / expected_hatch_bin_all);
    else
        keep_fraction = 0;
    end

    expected_hatch_all = expected_hatch_inc + keep_fraction * expected_hatch_bin_all;
    expected_worms_next = expected_survivors + expected_hatch_all;

    if total_eggs_bin > 0
        if deterministic_mode
            egg_bins_bin = keep_fraction * egg_bins_bin;
        else
            egg_bins_bin = binornd(round(egg_bins_bin), keep_fraction);
        end
    end

    % move kept bin eggs to incubator at reset
    egg_bins_inc = egg_bins_inc + egg_bins_bin;
    egg_bins_bin(:) = 0;

    if enable_debug
        fprintf('t=%6.1f d | horizon=%d cycles | keep=%.4f | worms now=%d | expected survivors=%.2f | expected hatch (controlled)=%.2f | expected worms next=%.2f\n', ...
            tk, control_horizon_cycles, keep_fraction, total_worms, expected_survivors, expected_hatch_all, expected_worms_next);
        fprintf('  Hatch schedule within horizon (time from now, absolute time, expected hatch count):\n');
        for i = 1:nEggBins
            age = (i-1) * dt;
            time_to_hatch = a_b_T - age;
            if time_to_hatch >= 0 && time_to_hatch <= horizon_len
                expected_hatch_i = egg_bins_inc(i) * exp(-mu_egg * time_to_hatch) * hatch_surv;
                if expected_hatch_i > 0
                    fprintf('    +%5.1f d (t=%.1f): %.2f eggs -> juveniles\n', ...
                        time_to_hatch, tk + time_to_hatch, expected_hatch_i);
                end
            end
        end
    end
end

keep_eggs = sum(egg_bins_inc);
keep_worms = numel(L_adult) + numel(L_juv);
end

function [q, h] = deb_aging_step(L, q, h, dLdt, f, L_m, v_T, h_a_T, s_G, dt)
% deb_aging_step: forward-Euler update for DEB aging acceleration and hazard; q and h, respectively.
% see Eqn. 6.2 of "Dynamic Energy Budget Theory for Metabolic Organisation" by S.A.L.M. Kooijman
L_safe = max(L, 1e-12);
r = dLdt ./ L_safe;
dq = (q .* s_G .* (L_safe.^3) ./ (L_m^3) + h_a_T) .* f .* (v_T ./ L_safe - r) - r .* q;
dh = q - r .* h;
q = max(q + dq .* dt, 0);
h = max(h + dh .* dt, 0);
end

function mask = deterministic_survival_mask(p_survive)
p_survive = p_survive(:);
n = numel(p_survive);
mask = false(n,1);
if n == 0
    return;
end
p_survive = min(max(p_survive, 0), 1);
n_keep = min(n, max(0, round(sum(p_survive))));
if n_keep <= 0
    return;
end
[~, idx] = sort(p_survive, 'descend');
mask(idx(1:n_keep)) = true;
end

function M = get_moisture(tk, moisture, moisture_t, moisture_series)
if ~isempty(moisture_t) && ~isempty(moisture_series)
    M = interp1(moisture_t, moisture_series, tk, 'linear', 'extrap');
else
    M = moisture;
end
M = min(max(M, 0), 0.99);
end

function [k_C_decay, k_N_loss, k_N_immob, frac_N_loss_assim, frac_N_immob_assim, extras] = ...
    eval_param_coeffs(coeffs, lb, ub, T_C, M, T_ref, M_ref, T_scale, M_scale, extra_lb, extra_ub, extra_defaults)

if nargin < 8 || isempty(T_scale); T_scale = 1.0; end
if nargin < 9 || isempty(M_scale); M_scale = 1.0; end
if nargin < 10 || isempty(extra_lb); extra_lb = [0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.00]; end
if nargin < 11 || isempty(extra_ub); extra_ub = [2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 0.08]; end
if nargin < 12 || isempty(extra_defaults)
    extra_defaults = struct('kC_sub_mult',1.0,'kC_faec_mult',1.0,'kNloss_sub_mult',1.0,'kNloss_faec_mult',1.0, ...
        'kNimmob_sub_mult',1.0,'kNimmob_faec_mult',1.0,'k_DM_other',0.0);
end

npar = size(coeffs, 1);
ncol = size(coeffs, 2);
b0 = coeffs(:,1);
bT = zeros(npar,1);
bM = zeros(npar,1);
bTM = zeros(npar,1);
if ncol >= 2; bT = coeffs(:,2); end
if ncol >= 3; bM = coeffs(:,3); end
if ncol >= 4; bTM = coeffs(:,4); end

xT = (T_C - T_ref) / max(T_scale, 1e-9);
xM = (M - M_ref) / max(M_scale, 1e-9);
theta = b0 + bT .* xT + bM .* xM + bTM .* (xT * xM);
s = 1 ./ (1 + exp(-theta));

lb_main_default = [1e-4; 1e-5; 1e-5; 0.00; 0.00];
ub_main_default = [0.20; 0.10; 0.10; 0.50; 0.90];
lb_main = lb(:);
ub_main = ub(:);
if numel(lb_main) < 5
    lb_main = [lb_main; lb_main_default(numel(lb_main)+1:5)];
else
    lb_main = lb_main(1:5);
end
if numel(ub_main) < 5
    ub_main = [ub_main; ub_main_default(numel(ub_main)+1:5)];
else
    ub_main = ub_main(1:5);
end

s_main = zeros(5,1);
s_main(1:min(5,npar)) = s(1:min(5,npar));
p_main = lb_main + s_main .* (ub_main - lb_main);
k_C_decay = p_main(1);
k_N_loss  = p_main(2);
k_N_immob = p_main(3);
frac_N_loss_assim  = p_main(4);
frac_N_immob_assim = p_main(5);

extras = extra_defaults;
s_extra = zeros(7,1);
if npar > 5
    n_extra = min(7, npar - 5);
    s_extra(1:n_extra) = s(6:(5+n_extra));
end
el = extra_lb(:);
eu = extra_ub(:);
if numel(el) < 7; el = [el; 0.40 * ones(7 - numel(el), 1)]; else; el = el(1:7); end
if numel(eu) < 7; eu = [eu; 2.50 * ones(7 - numel(eu), 1)]; else; eu = eu(1:7); end
eu(7) = max(eu(7), 0.0);
el(7) = max(el(7), 0.0);
p_extra = el + s_extra .* (eu - el);

extras.kC_sub_mult = p_extra(1);
extras.kC_faec_mult = p_extra(2);
extras.kNloss_sub_mult = p_extra(3);
extras.kNloss_faec_mult = p_extra(4);
extras.kNimmob_sub_mult = p_extra(5);
extras.kNimmob_faec_mult = p_extra(6);
extras.k_DM_other = p_extra(7);
end


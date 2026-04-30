function [prdData, info] = predict_Eisenia_fetida(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

  for i = 1:9
    ftW = eval(['ftW',num2str(i)]);
    if ftW > 1 
      info = 0; prdData = []; return
    end
  end 
% %   
  for i = 1:37
    ftW = eval(['ftW2',num2str(i)]);
    if ftW > 1 
      info = 0; prdData = []; return
    end
  end   

  % compute temperature correction factors
  TC_ab = tempcorr(temp.ab, T_ref, T_A);
  TC_tp = tempcorr(temp.tp, T_ref, T_A);
  TC_am = tempcorr(temp.am, T_ref, T_A);
  TC_Ri = tempcorr(temp.Ri, T_ref, T_A);
  TC_tW1 = tempcorr(temp.tW1, T_ref, T_A);
  TC_W11= tempcorr(temp.tW11, T_ref, T_A);

  % zero-variate data

  % life cycle
  pars_tp = [g k l_T v_Hb v_Hp];
  [t_p, t_b, l_p, l_b, info] = get_tp(pars_tp, f);
  
  % initial
%   pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
%   U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve

%   % hatch   
%   [U_H aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
%   a_h = aUL(2,1); aT_h = a_h/ TC_ah; % d, age at hatch at f and T

  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  Lw_b = L_b/ del_M;                % cm, physical length at birth at f
  Ww_b = L_b^3 * (1 + f * w);       % g, wet weight at birth at f 
  aT_b = t_b/ k_M/ TC_ab;           % d, age at birth at f and T

  % puberty 
  L_p = L_m * l_p;                  % cm, structural length at puberty at f
  Lw_p = L_p/ del_M;                % cm, physical length at puberty at f
  Ww_p = L_p^3 *(1 + f * w);        % g, wet weight at puberty 
  aT_p = (t_p - t_b)/ k_M/ TC_tp;           % d, age since birth at puberty at f and T

  % ultimate
  l_i = f - l_T;                    % -, scaled ultimate length
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Lw_i = L_i/ del_M;                % cm, ultimate physical length at f
  Ww_i = L_i^3 * (1 + f * w);       % g, ultimate wet weight 
 
  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp]; % compose parameter vector at T
  RT_i = TC_Ri * reprod_rate(L_i, f, pars_R);             % #/d, ultimate reproduction rate at T

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC_am;               % d, mean life span at T
  
  % pack to output
%   prdData.ah = aT_h;
  prdData.ab = aT_b;
  prdData.tp = aT_p;
  prdData.am = aT_m;
%   prdData.Lh = Lw_h;
  prdData.Lb = Lw_b;
  prdData.Lp = Lw_p;
  prdData.Li = Lw_i;
%   prdData.Wwh = Ww_h;
  prdData.Wwb = Ww_b;
  prdData.Wwp = Ww_p;
  prdData.Wwi = Ww_i;
  prdData.Ri  = RT_i;
  
  % uni-variate data : time-wet weight of 46 individual worms from 2
  % separate experiments
  
  L_i = L_m * (f - l_T);
  ir_B = 3/ k_M + 3 * f * L_m/ v; rT_B = TC_tW1/ ir_B;  % d, 1/von Bert growth rate

  % below we make a loop to compute the expected growth for all of the
  % invidividuals. This shorten the code cosiderably.
  % there are 9 individual worms in the first experiment and 37 worms in
  % the second experiment.
  
 % experiment 1: 
  for i = 1:9
    tW = ['tW',num2str(i)];
    prdData.(tW) = eval(['growth(tW',num2str(i),'(:,1),[ftW',num2str(i),'; W0.tW',num2str(i),'; w; L_i; rT_B])'] ) ; 
  end

% experiment 2: 
  for i = 1:37
    tW = ['tW2',num2str(i)];
    prdData.(tW) = eval(['growth(tW2',num2str(i),'(:,1),[ftW2',num2str(i),'; W0.tW2',num2str(i),'; w; L_i; rT_B])'] ) ;
  end
 
  
  % time-weight 
  f=ftW10; ir_B = 3/ k_M + 3 * f * L_m/ v;                          % d, 1/von Bertalanffy growth rate
  L_i = L_m * (f - l_T); L_b = L_m * get_lb([g k v_Hb], f);         % cm, structural length
  L = L_i - (L_i - L_b) * exp( - rT_B * tW10(:,1));                 % cm, structural length at time
  EWw10 = L.^3 * (1 + f * w);                                       % g, wet weight  
  % pack to output
  prdData.tW10 = EWw10;
  
  f=ftW11; ir_B = 3/ k_M + 3 * f * L_m/ v;  rT_B = TC_W11/ ir_B;    % d, 1/von Bertalanffy growth rate
  L_i = L_m * (f - l_T); L_b = L_m * get_lb([g k v_Hb], f);         % cm, structural length
  L = L_i - (L_i - L_b) * exp( - rT_B * tW11(:,1));                 % cm, structural length at time
  EWw11 = L.^3 * (1 + f * w);                                       % g, wet weight  
  % pack to output
  prdData.tW11 = EWw11;
  
  f=ftW12; ir_B = 3/ k_M + 3 * f * L_m/ v;  rT_B = TC_W11/ ir_B;    % d, 1/von Bertalanffy growth rate
  L_i = L_m * (f - l_T); L_b = L_m * get_lb([g k v_Hb], f);         % cm, structural length
  L = L_i - (L_i - L_b) * exp( - rT_B * tW12(:,1));                 % cm, structural length at time
  EWw12 = L.^3 * (1 + f * w);                                       % g, wet weight  
  % pack to output
  prdData.tW12 = EWw12;
  
  f=ftW13; ir_B = 3/ k_M + 3 * f * L_m/ v;  rT_B = TC_W11/ ir_B;    % d, 1/von Bertalanffy growth rate
  L_i = L_m * (f - l_T); L_b = L_m * get_lb([g k v_Hb], f);         % cm, structural length
  L = L_i - (L_i - L_b) * exp( - rT_B * tW13(:,1));                 % cm, structural length at time
  EWw13 = L.^3 * (1 + f * w);                                       % g, wet weight  
  % pack to output
  prdData.tW13 = EWw13;


%% subfunction for growth:

function EW = growth(tW, par)

f   = par(1);  % -, scaled func. response
W0  = par(2);  % g, initial observed wet weight
w   = par(3);  % -, contribution of reserve to weight
L_i = par(4);  % cm, ultimate length at f
r_B = par(5);  % r_B, von Bert growth rate at f and T

L_0 = (W0 / (1 + f * w)).^(1/ 3);
L = L_i - (L_i - L_0) * exp( - r_B * tW(:,1)); % cm, structural length
EW = L.^3 * (1 + w * f); % g, wet mass

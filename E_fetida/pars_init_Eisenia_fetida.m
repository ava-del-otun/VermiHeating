function [par, metaPar, txtPar] = pars_init_Eisenia_fetida(metaData)

metaPar.model = 'std'; 

%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 

%% core primary parameters 
par.z = 0.3431;       free.z     = 1;   units.z = '-';            label.z = 'zoom factor'; 
par.F_m = 6.5;        free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;      free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;      free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.01477;      free.v     = 1;   units.v = 'cm/d';         label.v = 'energy conductance'; 
par.kap = 0.5756;     free.kap   = 0;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.475;    free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M = 4740.8319;  free.p_M   = 1;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;          free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.002;      free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 4181.7688;  free.E_G   = 1;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 5.229e+00; free.E_Hb  = 1;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_Hp = 1.055e+03; free.E_Hp  = 1;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
par.h_a = 1.629e-08;  free.h_a   = 1;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 0.0001;     free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

%% other parameters 
par.T_A = 8000;       free.T_A   = 0;   units.T_A = 'K';          label.T_A = 'Arrhenius temperature'; 
par.del_M = 0.050632;  free.del_M = 1;   units.del_M = '-';        label.del_M = 'shape coefficient'; 
par.f = 1;            free.f     = 0;   units.f = '-';            label.f = 'scaled functional response for 0-var data'; 
par.ftW1 = 0.2125;    free.ftW1  = 0;   units.ftW1 = '-';         label.ftW1 = 'ind1exp1'; 
par.ftW10 = 1;        free.ftW10 = 0;   units.ftW10 = '-';        label.ftW10 = 'Gunadi'; 
par.ftW11 = 1;        free.ftW11 = 0;   units.ftW11 = '-';        label.ftW11 = 'Neuhauser sludge'; 
par.ftW12 = 0.9565;   free.ftW12 = 0;   units.ftW12 = '-';        label.ftW12 = 'Neuhauser cow manure'; 
par.ftW13 = 0.9356;   free.ftW13 = 0;   units.ftW13 = '-';        label.ftW13 = 'Neuhauser horse manure'; 
par.ftW2 = 0.2813;    free.ftW2  = 0;   units.ftW2 = '-';         label.ftW2 = 'ind2exp1'; 
par.ftW21 = 0.5001;   free.ftW21 = 0;   units.ftW21 = '-';        label.ftW21 = 'ind1exp1'; 
par.ftW210 = 0.3921;  free.ftW210 = 0;   units.ftW210 = '-';       label.ftW210 = 'ind10exp1'; 
par.ftW211 = 0.3524;  free.ftW211 = 0;   units.ftW211 = '-';       label.ftW211 = 'ind11exp1'; 
par.ftW212 = 0.5243;  free.ftW212 = 0;   units.ftW212 = '-';       label.ftW212 = 'ind12exp1'; 
par.ftW213 = 0.3943;  free.ftW213 = 0;   units.ftW213 = '-';       label.ftW213 = 'ind13exp1'; 
par.ftW214 = 0.4113;  free.ftW214 = 0;   units.ftW214 = '-';       label.ftW214 = 'ind14exp1'; 
par.ftW215 = 0.3988;  free.ftW215 = 0;   units.ftW215 = '-';       label.ftW215 = 'ind15exp1'; 
par.ftW216 = 0.3473;  free.ftW216 = 0;   units.ftW216 = '-';       label.ftW216 = 'ind16exp1'; 
par.ftW217 = 0.3655;  free.ftW217 = 0;   units.ftW217 = '-';       label.ftW217 = 'ind17exp1'; 
par.ftW218 = 0.188;   free.ftW218 = 0;   units.ftW218 = '-';       label.ftW218 = 'ind18exp1'; 
par.ftW219 = 0.2634;  free.ftW219 = 0;   units.ftW219 = '-';       label.ftW219 = 'ind19exp1'; 
par.ftW22 = 0.3672;   free.ftW22 = 0;   units.ftW22 = '-';        label.ftW22 = 'ind2exp1'; 
par.ftW220 = 0.17;    free.ftW220 = 0;   units.ftW220 = '-';       label.ftW220 = 'ind20exp1'; 
par.ftW221 = 0.1722;  free.ftW221 = 0;   units.ftW221 = '-';       label.ftW221 = 'ind21exp1'; 
par.ftW222 = 0.1958;  free.ftW222 = 0;   units.ftW222 = '-';       label.ftW222 = 'ind22exp1'; 
par.ftW223 = 0.1896;  free.ftW223 = 0;   units.ftW223 = '-';       label.ftW223 = 'ind23exp1'; 
par.ftW224 = 0.3862;  free.ftW224 = 0;   units.ftW224 = '-';       label.ftW224 = 'ind24exp1'; 
par.ftW225 = 0.3988;  free.ftW225 = 0;   units.ftW225 = '-';       label.ftW225 = 'ind25exp1'; 
par.ftW226 = 0.2758;  free.ftW226 = 0;   units.ftW226 = '-';       label.ftW226 = 'ind26exp1'; 
par.ftW227 = 0.2753;  free.ftW227 = 0;   units.ftW227 = '-';       label.ftW227 = 'ind27exp1'; 
par.ftW228 = 0.2058;  free.ftW228 = 0;   units.ftW228 = '-';       label.ftW228 = 'ind28exp1'; 
par.ftW229 = 0.3017;  free.ftW229 = 0;   units.ftW229 = '-';       label.ftW229 = 'ind29exp1'; 
par.ftW23 = 0.4038;   free.ftW23 = 0;   units.ftW23 = '-';        label.ftW23 = 'ind3exp1'; 
par.ftW230 = 0.25;    free.ftW230 = 0;   units.ftW230 = '-';       label.ftW230 = 'ind30exp1'; 
par.ftW231 = 0.3026;  free.ftW231 = 0;   units.ftW231 = '-';       label.ftW231 = 'ind31exp1'; 
par.ftW232 = 0.3161;  free.ftW232 = 0;   units.ftW232 = '-';       label.ftW232 = 'ind32exp1'; 
par.ftW233 = 0.3476;  free.ftW233 = 0;   units.ftW233 = '-';       label.ftW233 = 'ind33exp1'; 
par.ftW234 = 0.3037;  free.ftW234 = 0;   units.ftW234 = '-';       label.ftW234 = 'ind34exp1'; 
par.ftW235 = 0.28309;  free.ftW235 = 0;   units.ftW235 = '-';       label.ftW235 = 'ind35exp1'; 
par.ftW236 = 0.2979;  free.ftW236 = 0;   units.ftW236 = '-';       label.ftW236 = 'ind36exp1'; 
par.ftW237 = 0.2947;  free.ftW237 = 0;   units.ftW237 = '-';       label.ftW237 = 'ind37exp1'; 
par.ftW24 = 0.2947;   free.ftW24 = 0;   units.ftW24 = '-';        label.ftW24 = 'ind4exp1'; 
par.ftW25 = 0.512;    free.ftW25 = 0;   units.ftW25 = '-';        label.ftW25 = 'ind5exp1'; 
par.ftW26 = 0.4676;   free.ftW26 = 0;   units.ftW26 = '-';        label.ftW26 = 'ind6exp1'; 
par.ftW27 = 0.3813;   free.ftW27 = 0;   units.ftW27 = '-';        label.ftW27 = 'ind7exp1'; 
par.ftW28 = 0.3517;   free.ftW28 = 0;   units.ftW28 = '-';        label.ftW28 = 'ind8exp1'; 
par.ftW29 = 0.3088;   free.ftW29 = 0;   units.ftW29 = '-';        label.ftW29 = 'ind9exp1'; 
par.ftW3 = 0.2074;    free.ftW3  = 0;   units.ftW3 = '-';         label.ftW3 = 'ind3exp1'; 
par.ftW4 = 0.31516;   free.ftW4  = 0;   units.ftW4 = '-';         label.ftW4 = 'ind4exp1'; 
par.ftW5 = 0.2408;    free.ftW5  = 0;   units.ftW5 = '-';         label.ftW5 = 'ind5exp1'; 
par.ftW6 = 0.4859;    free.ftW6  = 0;   units.ftW6 = '-';         label.ftW6 = 'ind6exp1'; 
par.ftW7 = 0.2991;    free.ftW7  = 0;   units.ftW7 = '-';         label.ftW7 = 'ind7exp1'; 
par.ftW8 = 0.361;     free.ftW8  = 0;   units.ftW8 = '-';         label.ftW8 = 'ind8exp1'; 
par.ftW9 = 0.3578;    free.ftW9  = 0;   units.ftW9 = '-';         label.ftW9 = 'ind9exp1'; 

%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class);

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 

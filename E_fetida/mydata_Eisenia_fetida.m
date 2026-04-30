function [data, auxData, metaData, txtData, weights] = mydata_Eisenia_fetida 

%% set metaData
metaData.phylum     = 'Annelida'; 
metaData.class      = 'Clitellata'; 
metaData.order      = 'Haplotaxida'; 
metaData.family     = 'Lumbricidae';
metaData.species    = 'Eisenia_fetida'; 
metaData.species_en = 'Redworm'; 

metaData.ecoCode.climate = {'C'};
metaData.ecoCode.ecozone = {'TP', 'TN', 'TA', 'TH'};
metaData.ecoCode.habitat = {'0iTh', '0iTf', '0iTi'};
metaData.ecoCode.embryo  = {'Tt'};
metaData.ecoCode.migrate = {};
metaData.ecoCode.food    = {'biD'};
metaData.ecoCode.gender  = {'Hh'};
metaData.ecoCode.reprod  = {'O'};

metaData.T_typical  = C2K(20); % K, body temp
metaData.data_0     = {'ab'; 'ap'; 'am'; 'Lb'; 'Lp'; 'Li'; 'Wwb'; 'Wwp'; 'Wwi'; 'Ri'}; 
metaData.data_1     = {'t-Ww'}; 
metaData.COMPLETE = 2.5; % using criteria of LikaKear2011

metaData.author   = {'Starrlight Augustine', 'Dave Spurgeon', 'Ilze Rasnaca'};    
metaData.date_subm = [2016 09 22];              
metaData.email    = {'starrlight.augustine@akvaplan.niva.no'};            
metaData.address  = {'Akvaplan-niva, Fram High North Research Centre for Climate and the Environment, 9296 Tromso, Norway'}; 

metaData.author_mod_1    = {'Andre Gergs'};            
metaData.date_mod_1 = [2017 05 25];               
metaData.email_mod_1     = {'Gergs@gaiac.rwth-aachen.de'};            
metaData.address_mod_1   = {'gaiac, 52072 Aachen, Germany'};

metaData.curator     = {'Starrlight Augustine'};
metaData.email_cur   = {'sta@akvaplan.niva.no'}; 
metaData.date_acc    = [2017 07 22]; 


%% set data
% zero-variate data

data.ab = 25;    units.ab = 'd';    label.ab = 'age at birth';             bibkey.ab = 'CEH2016';   
  temp.ab = C2K(20);  units.temp.ab = 'K'; label.temp.ab = 'temperature';
data.tp = 60;    units.tp = 'd';    label.tp = 'time since birth at puberty';           bibkey.tp = 'CEH2016';
  temp.tp = C2K(20);  units.temp.tp = 'K'; label.temp.tp = 'temperature';
data.am = 5*365;    units.am = 'd';    label.am = 'life span';                bibkey.am = 'CEH2016';   
  temp.am = C2K(20);  units.temp.am = 'K'; label.temp.am = 'temperature'; 

data.Lb  = 1;   units.Lb  = 'cm';  label.Lb  = 'length at birth';   bibkey.Lb  = 'CEH2016';  
data.Lp  = 5;   units.Lp  = 'cm';  label.Lp  = 'length at puberty'; bibkey.Lp  = 'CEH2016'; 
data.Li  = 8;   units.Li  = 'cm';  label.Li  = 'ultimate total length';   bibkey.Li  = 'CEH2016';

% data.Wwb = 10*1e-3;   units.Wwb = 'g';   label.Wwb = 'wet weight at birth';     bibkey.Wwb = 'CEH2016';
data.Wwb = 15*1e-3;   units.Wwb = 'g';   label.Wwb = 'wet weight at birth';     bibkey.Wwb = 'ReinVent1987';
% data.Wwp = 275*1e-3;   units.Wwp = 'g';   label.Wwp = 'wet weight at puberty';   bibkey.Wwp = 'CEH2016';
data.Wwp = 600*1e-3;   units.Wwp = 'g';   label.Wwp = 'wet weight at puberty';   bibkey.Wwp = 'ReinVil1990';
% data.Wwi = 1000*1e-3;   units.Wwi = 'g';   label.Wwi = 'ultimate wet weight';     bibkey.Wwi = 'CEH2016';
data.Wwi = 2400*1e-3;   units.Wwi = 'g';   label.Wwi = 'ultimate wet weight';     bibkey.Wwi = 'Neuh1980';

data.Ri  = 2.7 * 2.3/7;   units.Ri  = '#/d'; label.Ri  = 'maximum reprod rate';     bibkey.Ri  = 'CEH2016';
comment.Ri = '2.3 cocoons/worm/week and average of 2.7 offspring per cocoon are produced. 2.7 individuals per concoon is in line with Reinecke and Kriel 19981, Dominguez(2004) reports higher numbers of individuals per cocoon (2.5-3.8)';
temp.Ri = C2K(20);    units.temp.Ri = 'K'; label.temp.Ri = 'temperature';

% data.Ri  = 4.6;   units.Ri  = '#/d'; label.Ri  = 'maximum reprod rate';     bibkey.Ri  = 'Hartenstein';
% comment.Ri = 'reproductive rate estimated cocoon production/individual per 15 weeks and from regressions: coocoon weight vs. adult weight and hatchlings/coocoon weight';
% temp.Ri = C2K(20);    units.temp.Ri = 'K'; label.temp.Ri = 'temperature'; 
% uni-variate data

% t-W data - experiment 1 (01 sep till 18 dec 2015)
% Based on CEH labwork by Ilze Rasnaca
% first weight is measured 28 days after the start of the experiment (ie
% when hatchlings of unknown age were placed in individual jars)
tW = [ ...
0	0.0266	0.0358	0.0029	0.0186	0.0082	0.0096	0.0052	0.0124	0.0078
31	0.0672	0.0858	0.0312	0.057	0.0392	0.0586	0.0252	0.0555	0.054
56	0.13	0.1765	0.0659	0.2009	0.1684	0.2985	0.1512	0.284	0.2761
84	0.2631	0.356	0.1747	0.2868	0.2319	0.4427	0.276	0.383	0.3763
108	0.3184	0.4053	0.2693	0.4019	0.2977	0.5522	0.3556	0.384	0.3446 ];
t = tW(:,1); % d, time
W = tW(:,2:end); % g, wet weight (each column is an individual work
% Here we creat a loop to make the data, and associated auxData for the
% nine individual worms:
for i = 1:length(W)
tW = ['tW',num2str(i)];
data.(tW) = [t, W(:,i)]; % d, gram
units.(tW)   = {'d', 'g'};  label.(tW) = {'time', 'wet weight'};  
temp.(tW)    = C2K(20);  units.temp.(tW) = 'K'; label.temp.(tW) = 'temperature';
auxData.W0.(tW) = W(1,i) ; units.W0.(tW)= 'g'; label.W0.(tW) = 'initial weight';
bibkey.(tW) = 'CEH2016';
end

% t-W data - experiment 2 
% Based on CEH labwork by Ilze Rasnaca
% hatchlings were put in jars some 6 weeks earlier than the start of the
% experiment
tW = [ ...
0	0.035	0.03	0.053	0.012	0.016	0.024	0.021	0.053	0.037	0.048	0.031	0.087	0.025	0.042	0.05	0.076	0.044	0.154	0.097	0.082	0.116	0.125	0.097	0.066	0.047	0.111	0.078	0.126	0.075	0.099	0.079	0.081	0.067	0.052	0.051	0.081	0.108
28	0.183	0.144	0.172	0.162	0.163	0.145	0.143	0.221	0.178	0.185	0.161	0.311	0.163	0.252	0.23	0.238	0.173	0.222	0.26	0.208	0.199	0.227	0.231	0.209	0.208	0.249	0.202	0.246	0.194	0.242	0.234	0.266	0.262	0.151	0.204	0.234	0.274
58	0.365	0.308	0.372	0.385	0.4	0.484	0.422	0.331	0.304	0.26	0.237	0.554	0.315	0.377	0.384	0.38	0.375	0.302	0.335	0.22	0.257	0.303	0.303	0.345	0.351	0.333	0.289	0.293	0.4	0.3	0.353	0.358	0.333	0.312	0.324	0.33	0.38
87	0.52	0.375	0.451	0.439	0.418	0.325	0.236	0.39	0.31	0.484	0.399	0.553	0.38	0.395	0.408	0.404	0.363	0.332	0.336	0.245	0.285	0.3	0.233	0.471	0.43	0.392	0.371	0.321	0.321	0.35	0.36	0.371	0.406	0.351	0.351	0.3745	0.3745
];
t = tW(:,1); % d, time
W = tW(:,2:end); % g, wet weight (each column is an individual worm
% Here we creat a loop to make the data, and associated auxData for the
% 37 individual worms:
for i = 1:length(W)
tW = ['tW2',num2str(i)];
data.(tW) = [t, W(:,i)]; % d, gram
units.(tW)   = {'d', 'g'};  label.(tW) = {'time', 'wet weight'};  
temp.(tW)    = C2K(20);  units.temp.(tW) = 'K'; label.temp.(tW) = 'temperature';
auxData.W0.(tW) = W(1,i) ; units.W0.(tW)= 'g'; label.W0.(tW) = 'initial weight';
bibkey.(tW) = 'CEH2016';
end

data.tW10 = [ ... % time since birth (d), wet weight (mg)
0	0.0109
7	0.032
14	0.079
21	0.131
28	0.256
35	0.368
42	0.414
49	0.442
56	0.462
% 63	0.448
% 70	0.450
% 77	0.459
% 84	0.494
% 91	0.495
% 98	0.490
];
units.tW10   = {'d', 'g'};  label.tW10 = {'time since birth', 'wet weight'}; 
temp.tW10    = C2K(20);  units.temp.tW10 = 'K'; label.temp.tW10 = 'temperature';
bibkey.tW10 = 'Guna2002'; comment.tW10='food was added once initially in the experiment and most probably depleted within the first ~60 days'; 

data.tW11 = [ ... % time since birth (d), wet weight (mg)
7	0.004412
14	0.052940
21	0.1456
28	0.2956
35	0.4853
42	0.7059
49	0.9176
56	1.0990
63	1.2840
70	1.3940
77	1.5040
84	1.7070
91	1.7910
98	2.0430
105	1.9680
112	2.1130
119	2.2320
126	2.3430
133	2.3650
140	2.3740
147	2.3740
154	2.3650
161	2.36
168	2.2810
175	2.3120
182	2.1180
189	2.1260
196	2.2150
203	2.2240
210	2.1790
217	2.1970
224	2.2540
231	2.0380
238	2.21
245	2.2280
252	2.1620];
units.tW11   = {'d', 'g'};  label.tW11 = {'time since birth', 'wet weight'}; 
temp.tW11    = C2K(25);  units.temp.tW11 = 'K'; label.temp.tW11 = 'temperature';
bibkey.tW11 = 'Neuh1980'; comment.tW11='food source: sludge';

% 
data.tW12 = [ ... % time since birth (d), wet weight (mg)
7	0.00882
14	0.01324
21	0.08824
28	0.1588
35	0.3044
42	0.5029
49	0.675
56	0.8471
63	0.9618
70	1.094
77	1.204
84	1.288
91	1.376
98	1.513
105	1.571
112	1.597
119	1.654
126	1.76
133	1.782
140	1.809	
154	1.866
161	1.888	
175	1.919	
189	1.831
196	1.862
203	1.844
210	1.857
217	1.901
224	1.857
231	1.866
238	1.95
245	1.932
252	1.968 
];
units.tW12   = {'d', 'g'};  label.tW12 = {'time since birth', 'wet weight'}; 
temp.tW12    = C2K(25);  units.temp.tW12 = 'K'; label.temp.tW12 = 'temperature';
bibkey.tW12 = 'Neuh1980'; comment.tW12='food source: cow manure';

% 
data.tW13 = [ ... % time since birth (d), wet weight (mg)
7	0
14	0
21	0.07059
28	0.1853
35	0.3309
42	0.4103
49	0.5029
56	0.5956
63	0.6971
70	0.8294
77	0.9441
84	1.041
91	1.134
98	1.204
105	1.249
112	1.293
119	1.341
126	1.429
133	1.482
140	1.562
154	1.557
161	1.69	
175	1.707
182	1.729
189	1.791
196	1.826
203	1.937
210	1.932
217	1.871
224	1.853
231	1.849
238	1.968
245	1.932
252	1.95];
units.tW13   = {'d', 'g'};  label.tW13 = {'time since birth', 'wet weight'}; 
temp.tW13    = C2K(25);  units.temp.tW13 = 'K'; label.temp.tW13 = 'temperature';
bibkey.tW13 = 'Neuh1980'; comment.tW13='food source: horse manure';

%% set weights for all real data
weights = setweights(data, []);

% --------- Do not delete these user defined lines: ------------------------
% for i = 1:9
% tW = ['tW',num2str(i)];
% weights.(tW) = weights.(tW) * 0.01; % d, gram
% end
% for i = 1:37
% tW = ['tW2',num2str(i)];
% weights.(tW) = weights.(tW) * 0.01; % d, gram
% end

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);
weights.psd.p_M = 10 * weights.psd.p_M;

%% pack auxData and txtData for output
auxData.temp = temp;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
txtData.comment = comment;

%% Group plots
set1 = {'tW1', 'tW2', 'tW3', 'tW4'}; subtitle1 = {'Weight for individuals, exp 1'};
set2 = {'tW5', 'tW6', 'tW7', 'tW8', 'tW9'}; subtitle2 = {'Weight for individuals, exp 1'};
set3 = {'tW21', 'tW22', 'tW23', 'tW24', 'tW25'}; subtitle3 = {'Weight for individuals, exp 2'};
set4 = {'tW26', 'tW27', 'tW28', 'tW29', 'tW210'}; subtitle4 = {'Weight for individuals, exp 2'};
set5 = {'tW211', 'tW212', 'tW213', 'tW214', 'tW215'}; subtitle5 = {'Weight for individuals, exp 2'};
set6 = {'tW216', 'tW217', 'tW218', 'tW219', 'tW220'}; subtitle6 = {'Weight for individuals, exp 2'};
set7 = {'tW221', 'tW222', 'tW223', 'tW224', 'tW225'}; subtitle7 = {'Weight for individuals, exp 2'};
set8 = {'tW226', 'tW227', 'tW228', 'tW229', 'tW230'}; subtitle8 = {'Weight for individuals, exp 2'};
set9 = {'tW231', 'tW232', 'tW233', 'tW234', 'tW235'}; subtitle9 = {'Weight for individuals, exp 2'};
set10 = {'tW236', 'tW237'}; subtitle10 = {'Weight for individuals, exp 2'};
% set11 = {'tW11', 'tW12', 'tW13'}; subtitle11 = {'Weight for differnt food sources'};
metaData.grp.sets = {set1, set2 , set3, set4, set5, set6 , set7, set8, set9, set10};
metaData.grp.subtitle = {subtitle1, subtitle2, subtitle3, subtitle4, subtitle5, subtitle6, subtitle7, subtitle8, subtitle9, subtitle10};

%% Discussion points
D1 = 'We assume kap_R 0 0.95/2 for consistency with how all hermaphrodites are treated in the collection';
D2 = 'The mod_1 included more weight against time since birth data.';     
metaData.discussion = struct('D1', D1, 'D2', D2);

%% Facts
F1 = 'Eisenia fetida is hermaphroditic. However, two worms are still required for reproduction.';
metaData.bibkey.F1 = 'Wiki'; 
metaData.facts = struct('F1',F1);

%% Acknowledgment
metaData.acknowledgment = 'The creation of this entry was supported by the European Food Safety Authority (grant number OC/EFSA/SCER/2015/01)';

%% Links
metaData.links.id_CoL = '38W46'; % Cat of Life
metaData.links.id_ITIS = '976620'; % ITIS
metaData.links.id_EoL = '3126934'; % Ency of Life
metaData.links.id_Wiki = 'Eisenia_fetida'; % Wikipedia
metaData.links.id_ADW = ''; % ADW
metaData.links.id_Taxo = '140603'; % Taxonomicon
metaData.links.id_WoRMS = '1041027'; % WoRMS

%% References
bibkey = 'Wiki'; type = 'Misc'; bib = ...
'howpublished = {\url{https://en.wikipedia.org/wiki/Eisenia_fetida}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Kooy2010'; type = 'Book'; bib = [ ...  % used in setting of chemical parameters and pseudodata
'author = {Kooijman, S.A.L.M.}, ' ...
'year = {2010}, ' ...
'title  = {Dynamic Energy Budget theory for metabolic organisation}, ' ...
'publisher = {Cambridge Univ. Press, Cambridge}, ' ...
'pages = {Table 4.2 (page 150), 8.1 (page 300)}, ' ...
'howpublished = {\url{../../../bib/Kooy2010.html}}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'CEH2016'; type = 'Misc'; bib = ...
'howpublished = {Centre for Hydrology and Ecology (CEH), 2016, pers. comm. Dave Spurgeon/Ilze Rasnaca}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Guna2002'; type = 'Article'; bib = [ ... 
'author = {Bintoro Gunadi, Charles Blount and Clive A. Edwards}, ' ... 
'year = {2002}, ' ...
'title = {The growth and fecundity of \emph{Eisenia fetida} ({S}avigny) in cattle solids pre-composted for different periods}, ' ...
'journal = {Pedobiologia (}, ' ...
'volume = {46}, ' ...
'pages = {15-23}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Neuh1980'; type = 'Article'; bib = [ ... 
'author = {Edward F. Neuhauser, Roy Hartenstein, David L. Kaplan}, ' ... 
'year = {1980}, ' ...
'title = {Growth of the Earthworm \emph{Eisenia foetida} in Relation to Population Density and Food Rationing}, ' ...
'journal = {Oikos}, ' ...
'volume = {35}, ' ...
'pages = {93-98}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'ReinVil1990'; type = 'Article'; bib = [ ... 
'author = {A. J. Reinecke and S. A. Viljoen}, ' ... 
'year = {1990}, ' ...
'title = {The influence of feeding patterns on growth and reproduction of the vermicomposting earthworm \emph{Eisenia fetida} ({O}ligochaeta)}, ' ...
'journal = {Biol Fertil Soils}, ' ...
'volume = {10}, ' ...
'pages = {184-187}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'ReinVent1987'; type = 'Article'; bib = [ ... 
'author = {A.J. Reinecke and J.M. Venter}, ' ... 
'year = {1987}, ' ...
'title = {Moisture preferences, growth and reproduction of the compost worm \emph{Eisenia fetida} ({O}ligochaeta)}, ' ...
'journal = {Biol Fertil Soils}, ' ...
'doi = {}, ' ...
'volume = {3}, ' ...
'pages = {135-141}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
 

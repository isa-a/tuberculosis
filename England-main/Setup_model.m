clear all;

states1     = {'U','Lf_imp','Lf','Ls','Pf_imp','Pf','Ps','Irem','Irec','Tx','Rlo','Rhi','R'};
gps.born    = {'dom','migr_rect','migr_long'};

[i, s, d, lim] = get_addresses({states1, gps.born}, [], [], [], 0);
d = char(d);

s.migr       = [s.migr_rect, s.migr_long];
s.allI       = [s.Irec, s.Irem];
s.infectious = s.allI;

% Include the auxiliaries
names = {'inc','incsources','mort','nTPT'};
lgths = [    3,           8,     1,     1];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    i.aux.(names{ii}) = inds;
    lim = inds(end);
end
i.nx = lim;

% -------------------------------------------------------------------------
% --- Set up selectors and aggregators

% --- Incidence
tmp = zeros(3,i.nstates); 
tmp(1, s.allI) = 1;
tmp(2, intersect(s.allI, s.migr)) = 1;                                     % Incidence from migrants
tmp(3, s.Irec)  = 1;                                                       % Incidence from 'recent transmission', i.e. Lf
agg.inc = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.allI,:) = 1;
% Remove transitions due to changing migrant status
tmp(s.migr_rect, s.migr_long) = 0;
tmp(s.migr_long, s.migr_rect) = 0;
sel.inc = tmp - diag(diag(tmp));


% --- Sources of incidence 
tmp = zeros(2,i.nstates);                                                  % Rows: 1.Dom all 2.Migr all, 3.Dom RR, 4.Migr RR
tmp(1,intersect(s.allI, s.dom)) = 1;
tmp(2,intersect(s.allI, s.migr)) = 1;
agg.incsources = sparse(tmp);

% Untreated TB infection
tmp = zeros(i.nstates);
tmp(s.allI, [s.Lf_imp, s.Lf, s.Ls]) = 1;
sel.L2I = sparse(tmp - diag(diag(tmp)));

% Treated TB infection
tmp = zeros(i.nstates);
tmp(s.allI, [s.Pf_imp, s.Pf, s.Ps]) = 1;
sel.P2I = sparse(tmp - diag(diag(tmp)));

% Relapse, non-post-treatment
tmp = zeros(i.nstates);
tmp(s.allI, [s.R, s.Rhi]) = 1;
sel.R2I = sparse(tmp - diag(diag(tmp)));

% Relapse, post-treatment
tmp = zeros(i.nstates);
tmp(s.allI, s.Rlo) = 1;
sel.T2I = sparse(tmp - diag(diag(tmp)));


% --- People starting TPT
tmp = zeros(i.nstates);
tmp([s.Pf_imp, s.Pf, s.Ps],[s.Lf_imp, s.Lf, s.Ls]) = 1;
sel.nTPT = tmp - diag(diag(tmp));




% -- Natural history parameters -------------------------------------------
r.progression0  = 0.0826;
r.LTBI_stabil0  = 0.872;
r.reactivation0 = 0.0006;

r.Tx         = 2;
r.ltfu       = (149/4365) / (3873/4365) * r.Tx;   
r.muTx       = (240/4365) / (3873/4365) * r.Tx;                                        
r.self_cure0 = 1/6;
r.muTB0      = 1/6;
r.relapse    = [0.032 0.14 0.0015];
p.imm        = 0.8;                                                        % Reduced susceptibility conferred by previous infection

p.LTBIinmigr   = 0.09;                                                     %ltbi inmigr 6pct

% --- Interventions 
p.migrTPT      = 0;                                                        % Proportion of migrants initiated on TPT on entry
p.TPTeff       = 0.81;                                                     % Effectiveness of TPT
r.TPT          = [0 0 0];                                                  % Uptake of TPT amongst: 1.domestic, 2.recent migrants, 3.long-term migrants
r.TPT2020rec   = 0.004;
r.ACF          = [0 0 0];
% r.ACF2         = [0 0 0 0];

% p.migrect_popn = 0.168;                                                  % Proportion of population that's migrant recent
r.migr         = 1/10;                                                     % Average duration of stay for migrants

% years = [2019, 2020, 2021, 2022, 2023];
% vals = [806,    682, 773, 1294, 1316];
% vals = vals * 1e5;
% pop = 68.35e6;
% prop = vals / pop;

years    = [2019, 2020, 2021, 2022, 2023, 2024, 2025];
vals     = [286, 210, 406, 734, 904, 509, 474]*1e3;
pop      = [66.6e6, 66.7e6, 67e6, 67.6e6, 68.2e6, 69.2e6, 69.4e6];

rin_vec  = vals./pop;

% -------------------------------------------------------------------------
% --- Name free parameters ------------------------------------------------

names = {'r_beta','p_betadec','r_gamma','r_migrout','rf_progr_migrVdom', 'p_migr_Lf', 'rf_muTB', 'rf_progression','rf_reactivation','rf_LTBIstabil','rf_selfcure', 'migrmix'};      
lgths =        [1,          1,        2,          1,                  1,           1,         1,                1,                1,              1,            1,         1];

lim = 0; xi = [];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    xi.(names{ii}) = inds;
    lim = inds(end);
end
xi.nx = lim;

% Set their boundaries
bds = [];
bds(xi.r_beta,:)            = [0 40];
bds(xi.p_betadec,:)         = [0 5];                                       % Percentage annual decline in beta
bds(xi.r_gamma,:)           = [0, 10; 1, 10];                              % Rate of treatment uptake, in 2015 and its relative value in 2020
bds(xi.r_migrout,:)         = [0 1];
bds(xi.rf_progr_migrVdom,:) = [1, 20];                                     % Relative rate of progression: migr_rect v domestic
bds(xi.p_migr_Lf,:)         = [0 1];                                       % Amongst newly-arriving migrants with LTBI, proportion having Lf   
bds(xi.rf_muTB,:)           = [0.5, 1.5];
bds(xi.rf_progression,:)    = [0.75, 1.25];
bds(xi.rf_reactivation,:)   = [0.75, 1.25];
bds(xi.rf_LTBIstabil,:)     = [0.75, 1.25];
bds(xi.rf_selfcure,:)       = [0.5, 1.5];
bds(xi.migrmix,:)           = [1 10];
prm.bounds                  = bds';


ref.i = i; ref.s = s; ref.xi = xi;
prm.p = p; prm.r = r; prm.agg = agg; prm.sel = sel;

% -------------------------------------------------------------------------
% --- Specify -------------------------------------------------------------

data.incd2010         = [15.6 16.6 17.6]; 
data.incd2020         = [7.3 8.36 9.3]; 
data.mort2020         = [0.28 0.3 0.32];                                   
data.p_migrTB         = [0.708 0.728 0.748];                               % Proportion contribution of migrants to overall incidence
data.p_migrpopn       = [0.138 0.168 0.198];                               % Proportion of population that is migrants
data.nTPT2019         = 1.3*[0.9 1 1.1];                                   % Number of TPT initiations in 2019, per 10^5 population
data.p_incd_recentinf = [20 25 30]/100;                                    
data.migr_incd_5yr    = [1, 2, 4];                                         % Cumulative 5-year incidence (PERCENT) of newly-arrived migrants

% data.p_LTBI_inmigr  = [0.1 0.15 0.2];                                     % Proportion of migrants with LTBI: from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8904125/
% data.p_LTBI_inmigr  = [0.05 0.06 0.07];                                   % Proportion of migrants with LTBI: from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8904125/
% data.p_vulnpopn     = [2 3.48 5.5]/100;                                   % Proportion of UK population being vulnerable
% data.p_vulnTB       = [6 9.05 12]/100;                                    % Proportion contribution to overall incidence
% data.incd_ch2020    = [0.5 1.3 2.5];                                      % incidence in children
% data.p_chpopn       = [0.1102 0.1902 0.2502];                             % proportion of country thats children
% data.p_adpopn       = [0.6154 0.8154 0.93];                               % proportion of country thats adults
% data.ch_notifs      = [101 141 181]/55.98e6*1e5;                          % notifications in the country


show = 0;
f1 = get_distribution_fns(data.incd2010,   'lognorm', show);
f2 = get_distribution_fns(data.incd2020,   'lognorm', show);
f3 = get_distribution_fns(data.mort2020,   'lognorm', show);
f4 = get_distribution_fns(data.p_migrTB,   'beta',    show);
f5 = get_distribution_fns(data.p_migrpopn, 'beta',    show);
f6 = get_distribution_fns(data.p_incd_recentinf, 'beta', show);
f7 = get_distribution_fns(data.migr_incd_5yr,    'lognorm', show);

% f5  = get_distribution_fns(data.p_LTBI_inmigr, 'beta',    show);
% f6  = get_distribution_fns(data.p_vulnpopn, 'beta',    show);
% f7  = get_distribution_fns(data.p_vulnTB,   'beta',    show);
% f8  = get_distribution_fns(data.incd_ch2020,'lognorm', show);
% f10  = get_distribution_fns(data.p_chpopn,  'beta',    show);
% f11  = get_distribution_fns(data.p_adpopn,  'beta',    show);
% f12  = get_distribution_fns(data.ch_notifs, 'lognorm', show);

lhd.fn = @(sim) f1(sim(1)) + f2(sim(2)) + f3(sim(3)) + f4(sim(4)) + f5(sim(5)) + f6(sim(6)) + f7(sim(7));

save Model_setup;
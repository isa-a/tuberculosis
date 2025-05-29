% UPDATE vuln data following most recent discussion with Sharon
% Tidy contact matrix and how it's implemented

clear all;

states1     = {'U'};
% states2     = {'Lf_imp','Lf','Ls','Pf','Ps','Irem','Irec','I2','Tx','Tx2','Rlo','Rhi','R'};
states2     = {'Lf_imp','Lf','Ls','Pf_imp','Pf','Ps','Irem','Irec','I2rem','I2rec','Tx','Tx2','Rlo','Rhi','R'};
gps.age     = {'ch','ad'};                                                  % added age here before other groups
% gps.born    = {'dom','migr_rect','migr_long','vuln'};
gps.born    = {'dom','migr_rect','migr_long'};
gps.strains = {'ds'};

[i, s, d, lim] = get_addresses({states1, gps.age, gps.born}, [], [], [], 0);
[i, s, d, lim] = get_addresses({states2, gps.age, gps.born, gps.strains}, i, s, d, lim);
d = char(d);

s.migr       = [s.migr_rect, s.migr_long];
s.allI       = [s.Irec, s.Irem, s.I2rec, s.I2rem];
% s.migrstates = [i.U.migr_rect, i.Lf.migr_rect, i.Ls.migr_rect, i.Pf.migr_rect, i.Ps.migr_rect];
% s.migrstates = intersect([s.U, s.Lf_imp, s.Lf, s.Ls, s.Pf, s.Ps], s.migr_rect);
s.infectious = s.allI;

% Include the auxiliaries
names = {'inc','incsources','mort','nTPT', 'ch_notifs'};
% lgths = [    5,          12,     1,     1,           1];
lgths = [    4,          12,     1,     1,           1];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    i.aux.(names{ii}) = inds;
    lim = inds(end);
end
i.nx = lim;


% -------------------------------------------------------------------------
% --- Set up selectors and aggregators

% --- Incidence
tmp = zeros(4,i.nstates); 
tmp(1,s.allI) = 1;
tmp(2,intersect(s.allI,s.migr)) = 1;
% tmp(3,intersect(s.allI,s.vuln)) = 1;
tmp(3,intersect(s.allI,s.ch))   = 1;
tmp(4,[s.Irec, s.I2rec])        = 1;
agg.inc = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.allI,:) = 1;
% Remove transitions due to changing migrant status
tmp(s.migr_rect, s.migr_long) = 0;
tmp(s.migr_long, s.migr_rect) = 0;
% Remove transitions due to changing drug resistance status

% Remove transitions due to change in vulnerability status
% tmp(s.dom, s.vuln) = 0;
% tmp(s.vuln, s.dom) = 0;
% Remove transitions due to change in age status
tmp(s.ch, s.ad) = 0;
tmp(s.ad, s.ch) = 0;
sel.inc = tmp - diag(diag(tmp));


% --- Sources of incidence 

% tmp = zeros(3,i.nstates);                                                  % Rows: 1.Dom all 2.Migr all, 3.Dom RR, 4.Migr RR
% tmp(1,intersect(s.allI, s.dom)) = 1;
% tmp(2,intersect(s.allI, s.migr)) = 1;
% tmp(3,intersect(intersect(s.allI, s.dom),s.rr)) = 1;
% tmp(4,intersect(intersect(s.allI, s.migr),s.rr)) = 1;
% agg.incsources = sparse(tmp);

% set1 = {s.dom, s.migr, s.vuln};
set1 = {s.dom, s.migr};
% set1 = {s.dom, s.migr};
set2 = {s.ds};

tmp  = zeros(length(set1)*length(set2),i.nstates);
row  = 1;
for is1 = 1:length(set1)
    for is2 = 1:length(set2)
        tmp(row, intersect(s.allI, intersect(set1{is1}, set2{is2}))) = 1;
        row = row+1;
    end
end
agg.incsources = sparse(tmp);


% --- Selectors for different origins 

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

% --- Notifications
tmp = zeros(i.nstates);
% tmp(intersect([s.Tx, s.Tx2],s.ch),:) = 1;
tmp(intersect(s.Tx,s.ch),:) = 1;
% Remove transitions due to changing migrant status
tmp(s.migr_rect, s.migr_long) = 0;
tmp(s.migr_long, s.migr_rect) = 0;
% Remove transitions due to changing drug resistance status

% Remove transitions due to change in vulnerability status
% tmp(s.dom, s.vuln) = 0;
% tmp(s.vuln, s.dom) = 0;
% Remove transitions due to change in age status
tmp(s.ch, s.ad) = 0;
tmp(s.ad, s.ch) = 0;
sel.ch_notifs = tmp - diag(diag(tmp));



% -- Natural history parameters -------------------------------------------
r.progression0  = 0.0826;
r.LTBI_stabil   = 0.872;
r.reactivation0 = 0.0006;

r.Tx            = 2;
r.ltfu          = 0.01;                                                     % loss to followup

% Second-line outcomes
% p.RRrec         = 1;
% r.RR_acqu       = 0.01;
% r.RR_acqu       = 0;
% r.Tx2           = 9/12;
% r.ltfu2         = r.Tx2*2;

r.self_cure    = 1/6;
r.relapse      = [0.032 0.14 0.0015];
%r.muTB         = 1/6;                                                      % TB related mortality
r.muTx         = -log(1 - 0.06)/0.5;                                           % mortality on tb treatment
p.imm          = 0.8;                                                      % Reduced susceptibility conferred by previous infection

p.ch_in_migr   = 0.07;                                                      % GUESS: check with countries

% --- Interventions 
p.migrTPT      = 0;                                                        % Proportion of migrants initiated on TPT on entry
% p.TPTeff       = [0.6 0.1];                                              % Effectiveness of TPT (DS and RR)
p.TPTeff       = 0.81;                                                     % Effectiveness of TPT
r.TPT          = [0 0 0 0];                                                % Uptake of TPT amongst: 1.domestic, 2.recent migrants, 3.long-term migrants, 4.vulnerables
r.TPT2020rec   = 0.004;
r.ACF          = [0 0 0 0];
r.ACF2         = [0 0 0 0];

p.migrect_popn = 0.168;                                                    % Proportion of population that's migrant recent
r.migr         = 1/10;                                                     % Average duration of stay for migrants

% p.LTBI_in_migrad = 0.17;
% p.LTBI_in_migrch = 0.03;

% years = [2019, 2020, 2021, 2022, 2023];
% vals = [806,    682, 773, 1294, 1316];
% vals = vals * 1e5;
% pop = 68.35e6;
% prop = vals / pop;

years    = [2019, 2020, 2021, 2022, 2023];
vals     = [806,    682, 773, 1294, 1316]*1e3;
pop      = [66.6e6, 66.7e6, 67e6, 67.6e6, 68.2e6];
rin_vec  = vals./pop;
% migrvec  = rin_vals/rin_vals(1);

% migrvec = [1.1792, 0.9978, 1.1309, 1.8932, 1.9254];
% migrvec2 = [1.2102, 1.0225, 1.1537, 1.9142, 1.9296];

% -------------------------------------------------------------------------
% --- Name free parameters ------------------------------------------------

% names = {'beta','betadec','gamma', 'p_birth','p_kLf','r_TPT2020rec','p_relrate','r_migr'};      
% lgths =      [1,        1,      1,         1,      1,             1,          1,       1];
% names = {'beta','betadec','gamma','r_TPT2020rec','p_relrate','r_migr','p_LTBI_in_migr'};      
% lgths =      [1,        1,      2,             1,          1,       1,               1];
% names = {'beta','betadec','gamma','p_relrate_gamma_chvad','p_relrate','r_migr','p_LTBI_in_migr','r_vuln','relbeta_vuln', 'ageing', 'ch_mort', 'p_relrate_factor'};      
% lgths =      [1,      1,      2,                  1,          2,            1,               1,       1,             1,       1,         1,                 1];
% names = {'beta','betadec','gamma','p_relrate_gamma_chvad','p_relrate','r_migr','p_LTBI_in_migr', 'ageing', 'ch_mort', 'p_relrate_factor'};      
% lgths =      [1,      1,      2,                  1,          2,            1,               1,        1,         1,                 1];
% names = {'beta','betadec','gamma','p_relrate_gamma_chvad','p_relrate','p_LTBI_in_migr', 'ageing', 'ch_mort', 'p_relrate_factor'};      
% lgths =      [1,      1,      2,                  1,          2,                  1,        1,         1,                 1];

% names = {'beta','betadec','gamma','p_relrate_gamma_chvad','p_LTBI_in_migrad','p_relLTBI_inmigr_advch','r_vuln_sc','relbeta_vuln', 'p_relrate', 'r_ageing_sc', 'p_relrate_factor'};      
% lgths =      [1,        1,      2,                      1,                 1,                       1,       1,             1,           2,        1,                  1];

names = {'beta','betadec','gamma','p_relrate_gamma_chvad','r_migrout','p_LTBI_in_migrad','p_relLTBI_inmigr_advch','p_relrate', 'r_ageing_sc', 'p_relrate_factor', 'muTB', 'migrmix'};      
lgths =      [1,        1,      2,                      1,          1,                 1,                       1,          1,             1,                  1,    1,           1];


lim = 0; xi = [];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    xi.(names{ii}) = inds;
    lim = inds(end);
end
xi.nx = lim;

% Set their boundaries
bds = [];
bds(xi.beta,:)                   = [0 40];
bds(xi.betadec,:)                = [0 0.15];
bds(xi.gamma,:)                  = [0, 1; 1e-4, 10];                       % Rate of treatment uptake, in 2015 and 2020
bds(xi.p_relrate_gamma_chvad,:)  = [0 2];                                  % Relative rate of treatment uptake, children vs adults
bds(xi.r_migrout,:)              = [0 1];
bds(xi.p_LTBI_in_migrad,:)       = [0.1 0.6];                              % Prevalence of LTBI in adult migrants, used to construct migrant 'birth' terms
bds(xi.p_relLTBI_inmigr_advch,:) = [3 7];                                  % Relative prev of LTBI in migrants, ad vs ch: estimated using simple ARTI calculation
% bds(xi.r_vuln_sc,:)              = [0 1];                                  % Rate of acquisition of vulnerability
% bds(xi.relbeta_vuln,:)           = [0.1 20];                               % Relative beta, vulnerable vs domestic - OVERWRITTEN
% bds(xi.p_relrate,:)              = repmat([1 20],2,1);                     % Relative rate of progression vs domestic, migr_recent and vulnerable 
bds(xi.p_relrate,:)              = [1 20];                                 % Relative rate of progression vs domestic, migr_recent and vulnerable 
bds(xi.r_ageing_sc,:)            = [0 1];                                  % Rate of ageing
bds(xi.p_relrate_factor,:)       = [1, 10];                                % Rel.rate of progression for children vs adults
%bds(xi.mort_factor,:)       = [0, 1];                               
bds(xi.muTB,:)       = [0, 1];
bds(xi.migrmix,:)       = [1 10];


scaling = ones(1,xi.nx);
% scaling(xi.r_vuln_sc)   = 500;
%scaling(xi.r_ageing_sc) = 10;
%scaling(xi.r_migrout) = 10;
prm.scaling             = scaling;

%bds(xi.r_migr,:)           = [0 1];
% bds(xi.r_vuln,:)                 = [0 2];
% bds(xi.ageing,:)                 = [0.02 0.3];
%bds(xi.ch_mort,:)          = [0, 0.01];
prm.bounds = bds';

ref.i = i; ref.s = s; ref.xi = xi;
prm.p = p; prm.r = r; prm.agg = agg; prm.sel = sel;

% prm.contmat_born = [1, 0.5, 0.2; 0.5, 1, 0.2; 0.2 0.2 1];
prm.contmat_born = [1, 0.5, 0.2; 0.5, 0.5, 0.2; 0.2, 0.2, 1];
prm.contmat_born(end,:) = [];
prm.contmat_born(:,end) = [];

prm.contmat_age  = [1 0.3; 0.3 1];

% prm.contmat      = zeros(6, 6);
prm.contmat      = zeros(4, 4);
% go through each element
for age_row = 1:2                                                           % rows in age
    for age_col = 1:2                                                       % cols in age
        for born_row = 1:2 %3                                                  % rows in born
            for born_col = 1:2 %3                                              % cols in born
                % calc position in combined matrix
                row = (born_row-1)*2 + age_row;                             % correctly scale current rows into new matrix
                col = (born_col-1)*2 + age_col;                             % correctly scale current cols into new matrix
                % multiply
                prm.contmat(row, col) = prm.contmat_born(born_row, born_col) * prm.contmat_age(age_row, age_col);
            end
        end
    end
end




% -------------------------------------------------------------------------
% --- Specify --------------------------------------------------------

% data.incd2010   = [14.1 14.6 15.1];
% data.incd2010       = [11 12 14];                                          % With broader uncertainty intervals
% data.incd2020       = [6.8  7.9  9.2];                                             
% data.incdRR2020     = [0.07 0.14 0.21];                                    % Incidence of RR-TB
% data.mort           = [0.26 0.36 0.47];                                    % TB mortality, 2020
% data.p_migrTB       = [0.5  0.6  0.7];                                     % Proportion contribution of migrants to overall incidence
% data.p_migrpopn     = [0.337 0.437 0.537];                                 % Proportion of population that is migrants
% data.p_LTBI         = [0.18 0.22 0.28];                                    % Proportion of migrants with LTBI: from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8904125/
% data.p_vulnpopn     = [5 10 15]/100;                                       % Proportion of UK population being vulnerable
% data.p_vulnTB       = [13 18 23]/100;                                      % Proportion contribution to overall incidence
% data.nTPT2019       = 1.3*[0.9 1 1.1];                                     % Number of TPT initiations in 2019, per 10^5 population
% data.incd_ch2020    = [0.5 1.3 2.5];                                       % incidence in children
% data.p_chpopn       = [0.286 0.292 0.3];                                   % proportion of country thats children
% data.p_adpopn       = [0.688 0.708 0.728];                                 % proportion of country thats adults
% data.ch_notifs      = [310 360 420]/4.5e6*1e5;                             % notifications in the country   

data.incd2010       = [13.6 14.6 15.6];                                    % With broader uncertainty intervals
data.incd2020       = [6 7 8];                                             
%data.incdRR2020     = [0.11 0.15 0.18];                                    % Incidence of RR-TB
data.mort           = [0.28 0.3 0.32];                                     % TB mortality, 2020
data.p_migrTB       = [0.708 0.728 0.748];                                 % Proportion contribution of migrants to overall incidence
data.p_migrpopn     = [0.138 0.168 0.198];                                 % Proportion of population that is migrants
data.p_LTBI_inmigr  = [0.15 0.2 0.25];                                     % Proportion of migrants with LTBI: from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8904125/
data.p_vulnpopn     = [2 3.48 5.5]/100;                                         % Proportion of UK population being vulnerable
data.p_vulnTB       = [6 9.05 12]/100;                                        % Proportion contribution to overall incidence
data.nTPT2019       = 1.3*[0.9 1 1.1];                                     % Number of TPT initiations in 2019, per 10^5 population
data.incd_ch2020    = [0.5 1.3 2.5];                                       % incidence in children
data.p_chpopn       = [0.1102 0.1902 0.2502];                              % proportion of country thats children
data.p_adpopn       = [0.6154 0.8154 0.93];                                % proportion of country thats adults
data.ch_notifs      = [101 141 181]/55.98e6*1e5;                           % notifications in the country
data.p_incd_recentinf = [20 25 30]/100;                                    

show = 0;
f1a = get_distribution_fns(data.incd2010,   'lognorm', show);
f1b = get_distribution_fns(data.incd2020,   'lognorm', show);
%f1c = get_distribution_fns(data.incdRR2020, 'lognorm', show);
f2  = get_distribution_fns(data.mort,       'lognorm', show);
f3  = get_distribution_fns(data.p_migrTB,   'beta',    show);
f4  = get_distribution_fns(data.p_migrpopn, 'beta',    show);
f5  = get_distribution_fns(data.p_LTBI_inmigr, 'beta',    show);
f6  = get_distribution_fns(data.p_vulnpopn, 'beta',    show);
f7  = get_distribution_fns(data.p_vulnTB,   'beta',    show);
f8  = get_distribution_fns(data.incd_ch2020,'lognorm', show);
f10  = get_distribution_fns(data.p_chpopn,  'beta',    show);
f11  = get_distribution_fns(data.p_adpopn,  'beta',    show);
f12  = get_distribution_fns(data.ch_notifs, 'lognorm', show);
f13  = get_distribution_fns(data.p_incd_recentinf, 'beta', show);

% lhd.fn = @(incd, mort, p_migrTB, p_migrpopn, p_LTBI) f1(incd) + f2(mort) + f3(p_migrTB) + f4(p_migrpopn) + f5(p_LTBI);
% lhd.fn = @(incd2010, incd2020, mort, p_migrTB, p_migrpopn, p_LTBI, nTPT2019) f1a(incd2010) + f1b(incd2020) + f2(mort) + f3(p_migrTB) + f4(p_migrpopn) + f5(p_LTBI) + f6(nTPT2019);

% lhd.fn = @(incd2010, incd2020, mort, p_migrTB, p_LTBI_inmigr, p_vulnpopn, p_vulnTB, incd_ch2020, p_chpopn, p_adpopn, ch_notifs) f1a(incd2010) + f1b(incd2020) + ...
%                                                                                                      f2(mort) + f3(p_migrTB) + f5(p_LTBI_inmigr) + ...
%                                                                                                      f6(p_vulnpopn) + f7(p_vulnTB) + ...
%                                                                                                      f8(incd_ch2020) + ...
%                                                                                                      f10(p_chpopn) + f11(p_adpopn) + f12(ch_notifs);

lhd.fn = @(incd2010, incd2020, mort, p_migrTB, p_migrpopn, p_LTBI_inmigr, p_chpopn, ch_notifs, p_incd_recentinf) f1a(incd2010) + f1b(incd2020) + ...
                                                                                                     f2(mort) + f3(p_migrTB) + f4(p_migrpopn) + f5(p_LTBI_inmigr) + ...                                                                                                    
                                                                                                     f10(p_chpopn) + f12(ch_notifs) + f13(p_incd_recentinf);



save Model_setup;
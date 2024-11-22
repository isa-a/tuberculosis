% UPDATE vuln data following most recent discussion with Sharon
% Tidy contact matrix and how it's implemented

clear all;

states1     = {'U'};
states2     = {'Lf','Ls','Pf','Ps','I','I2','Tx','Tx2','Rlo','Rhi','R'};
gps.age     = {'ch','ad','el'};                                                  % added age here before other groups
gps.born    = {'dom','migr_rect','migr_long','vuln'};
gps.strains = {'ds'};

[i, s, d, lim] = get_addresses({states1, gps.age, gps.born}, [], [], [], 0);
[i, s, d, lim] = get_addresses({states2, gps.age, gps.born, gps.strains}, i, s, d, lim);
d = char(d);

s.migr       = [s.migr_rect, s.migr_long];
s.allI       = [s.I, s.I2];
% s.migrstates = [i.U.migr_rect, i.Lf.migr_rect, i.Ls.migr_rect, i.Pf.migr_rect, i.Ps.migr_rect];
s.migrstates = intersect([s.U, s.Lf, s.Ls, s.Pf, s.Ps],s.migr_rect);
s.infectious = [s.allI, (s.Tx)];

% Include the auxiliaries
names = {'inc','incsources','mort','nTPT', 'ch_notifs'};
lgths = [    5,          12,     1,     1,           1];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    i.aux.(names{ii}) = inds;
    lim = inds(end);
end
i.nx = lim;


% -------------------------------------------------------------------------
% --- Set up selectors and aggregators

% --- Incidence
tmp = zeros(5,i.nstates); 
tmp(1,s.allI) = 1;
tmp(2,intersect(s.allI,s.migr)) = 1;
tmp(3,intersect(s.allI,s.vuln)) = 1;
tmp(4,intersect(s.allI,s.ch))   = 1;
tmp(5,intersect(s.allI,s.el))   = 1;
agg.inc = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.allI,:) = 1;
% Remove transitions due to changing migrant status
tmp(s.migr_rect, s.migr_long) = 0;
tmp(s.migr_long, s.migr_rect) = 0;
% Remove transitions due to changing drug resistance status

% Remove transitions due to change in vulnerability status
tmp(s.dom, s.vuln) = 0;
tmp(s.vuln, s.dom) = 0;
% Remove transitions due to change in age status
tmp(s.ch, s.ad) = 0;
tmp(s.ad, s.ch) = 0;

tmp(s.ad, s.el) = 0;
tmp(s.el, s.ad) = 0;
sel.inc = tmp - diag(diag(tmp));


% --- Sources of incidence 

% tmp = zeros(3,i.nstates);                                                  % Rows: 1.Dom all 2.Migr all, 3.Dom RR, 4.Migr RR
% tmp(1,intersect(s.allI, s.dom)) = 1;
% tmp(2,intersect(s.allI, s.migr)) = 1;
% tmp(3,intersect(intersect(s.allI, s.dom),s.rr)) = 1;
% tmp(4,intersect(intersect(s.allI, s.migr),s.rr)) = 1;
% agg.incsources = sparse(tmp);

set1 = {s.dom, s.migr, s.vuln};
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
tmp(s.allI, [s.Lf, s.Ls]) = 1;
sel.L2I = sparse(tmp - diag(diag(tmp)));

% Treated TB infecction
tmp = zeros(i.nstates);
tmp(s.allI, [s.Pf, s.Ps]) = 1;
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
tmp([s.Pf, s.Ps],[s.Lf, s.Ls]) = 1;
sel.nTPT = tmp - diag(diag(tmp));

% --- Notifications
tmp = zeros(i.nstates);
tmp(intersect([s.Tx, s.Tx2],s.ch),:) = 1;
% Remove transitions due to changing migrant status
tmp(s.migr_rect, s.migr_long) = 0;
tmp(s.migr_long, s.migr_rect) = 0;
% Remove transitions due to changing drug resistance status

% Remove transitions due to change in vulnerability status
tmp(s.dom, s.vuln) = 0;
tmp(s.vuln, s.dom) = 0;
% Remove transitions due to change in age status
tmp(s.ch, s.ad) = 0;
tmp(s.ad, s.ch) = 0;

tmp(s.ad, s.el) = 0;
tmp(s.el, s.ad) = 0;
sel.ch_notifs = tmp - diag(diag(tmp));


% -- Natural history parameters -------------------------------------------
r.progression0  = 0.0826;
r.LTBI_stabil   = 0.872;
r.reactivation0 = 0.0006;

r.Tx            = 2;
p.RRrec         = 1;
r.RR_acqu       = 0;
r.Tx2           = 9/12;
r.ltfu          = 0.01;                                                     % loss to followup
r.ltfu2         = r.Tx2*2;

r.self_cure    = 1/6;
r.relapse      = [0.032 0.14 0.0015];
% r.relapse      = [0 0 0];
% r.mu           = 1/66;                                                   % natural mortality
r.muTB         = 1/6;                                                      % TB related mortality
p.imm          = 0.8;                                                      % Reduced susceptibility conferred by previous infection

p.ch_in_migr   = 0.0789;                                                      % DOUBLE check with country data

% --- Interventions 
p.migrTPT      = 0;                                                        % Proportion of migrants initiated on TPT on entry
% p.TPTeff       = [0.6 0.1];                                                % Effectiveness of TPT
p.TPTeff       = 0.6;                                                % Effectiveness of TPT
r.TPT          = [0 0 0 0];                                                % Uptake of TPT amongst: 1.domestic, 2.recent migrants, 3.long-term migrants
r.TPT2020rec   = 0.004;
r.ACF          = [0 0 0 0];
r.ACF2         = [0 0 0 0];

p.migrect_popn = 0.0419;
r.migr         = 0.0847;                                                   % https://doi.org/10.1007/s44197-022-00040-w

% p.LTBI_in_migrad = 0.17;
% p.LTBI_in_migrch = 0.03;


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

names = {'beta','betadec','gamma','p_relrate_gamma_chvad','p_LTBI_in_migrad','p_relLTBI_inmigr_advch','r_vuln_sc','relbeta_vuln', 'p_relrate', 'r_ageing_sc','p_relrate_factor', 'contmat_factor', 'ageing_ad2el'};      
lgths =      [1,        1,      2,                      1,                 1,                       1,       1,             1,           2,        1,            1,                1,                           1];

lim = 0; xi = [];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    xi.(names{ii}) = inds;
    lim = inds(end);
end

% Set their boundaries
bds = [];
bds(xi.beta,:)                   = [0 40];
bds(xi.betadec,:)                = [0 0.15];
bds(xi.gamma,:)                  = repmat([1e-4 10],2,1);
bds(xi.p_relrate_gamma_chvad,:)  = [0 1];
bds(xi.p_LTBI_in_migrad,:)       = [0.1 0.6];
bds(xi.p_relLTBI_inmigr_advch,:) = [3 7];                                % Estimated using simple ARTI calculation
bds(xi.p_relrate,:)              = repmat([1 20],2,1);
%bds(xi.r_migr,:)           = [0 1];
% bds(xi.r_vuln,:)                 = [0 2];
bds(xi.r_vuln_sc,:)              = [0 1];
bds(xi.relbeta_vuln,:)           = [0.1 20];
% bds(xi.ageing,:)                 = [0.02 0.3];
bds(xi.r_ageing_sc,:)            = [0 1];
%bds(xi.ch_mort,:)          = [0, 0.01];
bds(xi.p_relrate_factor,:) = [1, 10];
bds(xi.contmat_factor,:)    = [0, 1];
bds(xi.ageing_ad2el,:)    = [0.02, 3];
prm.bounds = bds';

ref.i = i; ref.s = s; ref.xi = xi;
prm.p = p; prm.r = r; prm.agg = agg; prm.sel = sel;

prm.contmat_born = [1, 0.5, 0.2; 0.5, 1, 0.2; 0.2 0.2 1];
prm.contmat_age  = [0.2830 0.2525 0.2; 0.0692 0.3953 0.2; 0.0692 0.3953 0.2];

prm.contmat      = zeros(9, 9);
% go through each element
for age_row = 1:3                                                           % rows in age
    for age_col = 1:3                                                      % cols in age
        for born_row = 1:3                                                  % rows in born
            for born_col = 1:3                                              % cols in born
                % calc position in combined matrix
                row = (born_row-1)*3 + age_row;                             % correctly scale current rows into new matrix
                col = (born_col-1)*3 + age_col;                             % correctly scale current cols into new matrix
                % multiply
                prm.contmat(row, col) = prm.contmat_born(born_row, born_col) * prm.contmat_age(age_row, age_col);
            end
        end
    end
end




% -------------------------------------------------------------------------
% --- Specify --------------------------------------------------------


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

data.incd2010       = [15 18 21];                                           % With broader uncertainty intervals
data.incd2020       = [10 12  14];                                             
%data.incdRR2020     = [0.11 0.15 0.18];                                    % Incidence of RR-TB
data.mort           = [0.8 1.1 1.5];                                     % TB mortality, 2020
data.p_migrTB       = [0.5  0.6  0.7];                                      % Proportion contribution of migrants to overall incidence
data.p_migrpopn     = [0.02 0.0419 0.07];                                  % Proportion of population that is migrants
data.p_LTBI_inmigr  = [0.16 0.20 0.24];                                     % Proportion of migrants with LTBI: from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8904125/
data.p_vulnpopn     = [8 10 12]/100;                                        % Proportion of UK population being vulnerable
data.p_vulnTB       = [5 7 9]/100;                                          % Proportion contribution to overall incidence
data.nTPT2019       = 1.3*[0.9 1 1.1];                                      % Number of TPT initiations in 2019, per 10^5 population
%data.incd_ch2020    = [0.05 0.11 0.16];                                    % incidence in children
data.propincd_ch    = [0.006 0.014 0.025];
data.p_chpopn       = [0.188 0.224 0.29];                                    % proportion of country thats children
data.p_adpopn       = [0.675 0.694 0.8275];                                  % proportion of country thats adults
data.p_elpopn       = [0.675 0.082 0.8275];                                  % proportion of country thats elderly
data.ch_notifs      = [410 450 490]/91567738*1e5;                             % notifications in the country  
% data.vuln_prev      = [10 12 14]/100;                                     
% data.vuln_relrisk   = [2 3 4];                                         


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
%f8  = get_distribution_fns(data.incd_ch2020,'lognorm', show);
f9  = get_distribution_fns(data.propincd_ch,'beta', show);
f10  = get_distribution_fns(data.p_chpopn,  'beta',    show);
f11  = get_distribution_fns(data.p_adpopn,  'beta',    show);
f12  = get_distribution_fns(data.p_elpopn,  'beta',    show);
f13  = get_distribution_fns(data.ch_notifs, 'lognorm', show);
% f13  = get_distribution_fns(data.vuln_prev, 'beta', show);
% f14  = get_distribution_fns(data.vuln_relrisk, 'lognorm', show);

% lhd.fn = @(incd, mort, p_migrTB, p_migrpopn, p_LTBI) f1(incd) + f2(mort) + f3(p_migrTB) + f4(p_migrpopn) + f5(p_LTBI);
% lhd.fn = @(incd2010, incd2020, mort, p_migrTB, p_migrpopn, p_LTBI, nTPT2019) f1a(incd2010) + f1b(incd2020) + f2(mort) + f3(p_migrTB) + f4(p_migrpopn) + f5(p_LTBI) + f6(nTPT2019);
lhd.fn = @(incd2010, incd2020, mort, p_migrTB, p_LTBI_inmigr, p_vulnpopn, p_vulnTB, p_chpopn, ch_notifs) f1a(incd2010) + f1b(incd2020) + ...
                                                                                                     f2(mort) + f3(p_migrTB) + f5(p_LTBI_inmigr) + ...
                                                                                                     f6(p_vulnpopn) + f7(p_vulnTB) + f10(p_chpopn) +f12(ch_notifs);

save Model_setup;


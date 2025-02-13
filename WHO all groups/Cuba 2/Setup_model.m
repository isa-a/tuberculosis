clear all;

states1     = {'U'};
states2     = {'Lf','Ls','Pf','Ps','I','I2','Tx','Tx2','Rlo','Rhi','R'};
gps.age     = {'ch','ad'};                                                  % added age here before other groups
gps.born    = {'dom','vuln'};
gps.strains = {'ds'};
gps.hiv     = {'neg','pos','art'};

[i, s, d, lim] = get_addresses({states1, gps.age, gps.born, gps.hiv}, [], [], [], 0);
[i, s, d, lim] = get_addresses({states2, gps.age, gps.born, gps.strains, gps.hiv}, i, s, d, lim);
d = char(d);

%s.migr       = [s.migr_rect, s.migr_long];
s.allI       = [s.I, s.I2];
% s.migrstates = [i.U.migr_rect, i.Lf.migr_rect, i.Ls.migr_rect, i.Pf.migr_rect, i.Ps.migr_rect];
%s.migrstates = intersect([s.U, s.Lf, s.Ls, s.Pf, s.Ps],s.migr_rect);
s.infectious = [s.allI, (s.Tx)];
s.prev       = [s.allI, (s.Tx)];

% Include the auxiliaries
names = {'inc','incsources','mort','nTPT', 'ch_notifs'};
lgths = [    6,          36,     1,     1,           1];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    i.aux.(names{ii}) = inds;
    lim = inds(end);
end
i.nx = lim;


% -------------------------------------------------------------------------
% --- Set up selectors and aggregators

% --- Incidence
tmp = zeros(6,i.nstates); 
tmp(1,s.allI) = 1;
tmp(2,intersect(s.allI,s.vuln)) = 1;
tmp(3,intersect(s.allI,s.ch))   = 1;
tmp(4,intersect(s.allI,s.neg))   = 1;
tmp(5,intersect(s.allI,s.pos))   = 1;
tmp(6,intersect(s.allI,s.art))   = 1;
agg.inc = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.allI,:) = 1;
% Remove transitions due to change in vulnerability status
tmp(s.dom, s.vuln) = 0;
tmp(s.vuln, s.dom) = 0;
% Remove transitions due to change in age status
tmp(s.ch, s.ad) = 0;
tmp(s.ad, s.ch) = 0;
% Remove transitions due to change in HIV status
tmp(s.pos, s.neg) = 0;
tmp(s.neg, s.pos) = 0;
tmp(s.art,s.pos) = 0; 
tmp(s.pos,s.art) = 0;
sel.inc = tmp - diag(diag(tmp));


% --- Sources of incidence 

set1 = {s.dom, s.vuln};
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
% Remove transitions due to change in vulnerability status
tmp(s.dom, s.vuln) = 0;
tmp(s.vuln, s.dom) = 0;
% Remove transitions due to change in age status
tmp(s.ch, s.ad) = 0;
tmp(s.ad, s.ch) = 0;
% Remove transitions due to change in HIV status
tmp(s.pos, s.neg) = 0;
tmp(s.neg, s.pos) = 0;
tmp(s.art,s.pos) = 0; 
tmp(s.pos,s.art) = 0;
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

%p.ch_in_migr   = 0.0789;                                                      % DOUBLE check with country data

% --- Interventions 
p.migrTPT      = 0;                                                        % Proportion of migrants initiated on TPT on entry
% p.TPTeff       = [0.6 0.1];                                                % Effectiveness of TPT
p.TPTeff       = 0.6;                                                % Effectiveness of TPT
r.TPT          = [0 0 0 0];                                                % Uptake of TPT amongst: 1.domestic, 2.recent migrants, 3.long-term migrants
r.TPT2020rec   = 0.004;
r.ACF          = [0 0 0 0];
r.ACF2         = [0 0 0 0];

% p.migrect_popn = 0.437;
% r.migr         = 0.0847;                                                   % https://doi.org/10.1007/s44197-022-00040-w

% p.LTBI_in_migrad = 0.17;
% p.LTBI_in_migrch = 0.03;


% -------------------------------------------------------------------------
% --- Name free parameters ------------------------------------------------

names = {'beta','betadec','gamma','p_relrate_gamma_chvad','r_vuln_sc','relbeta_vuln', 'p_relrate', 'r_ageing_sc','p_relrate_factor', 'contmat_factor', 'r_ART_init','r_HIV_mort','p_HIV_relrate'};      
lgths = [    1,        1,      2,                 1,        1,             1,           2,              1,                  1,                1,            1,              1,              1];

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
bds(xi.p_relrate,:)              = repmat([1 20],2,1);
bds(xi.r_vuln_sc,:)              = [0 1];
bds(xi.relbeta_vuln,:)           = [0.1 20];
bds(xi.r_ageing_sc,:)            = [0 1];
bds(xi.p_relrate_factor,:) = [1, 10];
bds(xi.contmat_factor,:)    = [0, 1];
bds(xi.r_ART_init,:)    = [0 10];
bds(xi.r_HIV_mort,:)    = [0 10];
bds(xi.p_HIV_relrate,:) = [1 100];
prm.bounds = bds';

ref.i = i; ref.s = s; ref.xi = xi;
prm.p = p; prm.r = r; prm.agg = agg; prm.sel = sel;

prm.contmat_born = [1, 0.5, 0.2; 0.5, 1, 0.2; 0.2 0.2 1];
prm.contmat_age  = [0.2830 0.2525; 0.0692 0.3953];

prm.contmat      = zeros(4, 4);
% go through each element
for age_row = 1:2                                                           % rows in age
    for age_col = 1:2                                                       % cols in age
        for born_row = 1:2                                                  % rows in born
            for born_col = 1:2                                              % cols in born
                % calc position in combined matrix
                row = (born_row-1)*2 + age_row;                             % correctly scale current rows into new matrix
                col = (born_col-1)*2 + age_col;                             % correctly scale current cols into new matrix
                % multiply
                prm.contmat(row, col) = prm.contmat_born(born_row, born_col) * prm.contmat_age(age_row, age_col);
            end
        end
    end
end



load data_Thembisa_AFR.mat;

% -------------------------------------------------------------------------
% --- Specify --------------------------------------------------------

data.incd2010       = [10 12 14];                                           % With broader uncertainty intervals
data.incd2020       = [6.8  7.9  9.2];                                             
data.mort           = [0.26 0.36 0.47];                                     % TB mortality, 2020
data.p_vulnpopn     = [8 10 12]/100;                                        % Proportion of UK population being vulnerable
data.p_vulnTB       = [5 7 9]/100;                                          % Proportion contribution to overall incidence
data.nTPT2019       = 1.3*[0.9 1 1.1];                                      % Number of TPT initiations in 2019, per 10^5 population
data.propincd_ch    = [0.006 0.014 0.025];
data.p_chpopn       = [0.198 0.2471 0.3];                                    % proportion of country thats children
data.p_adpopn       = [0.65 0.7529 0.85];                                  % proportion of country thats adults
data.ch_notifs      = [3 3.69 4.2]/4.576e6*1e5;                             % notifications in the country  
data.HIV_prev       = HIVprev_2019;
data.ART_covg       = ARTcovg_2019;
prm.ART_start       = ART_start;

ys1 = HIV_incd(:,2);
xs1 = 1:length(ys1);
xs2 = 1:length(ys1)+30;
ys2 = max(interp1(xs1,ys1,xs2,'linear','extrap'),0);
prm.rHIV = ys2;

show = 0;
f1a = get_distribution_fns(data.incd2010,   'lognorm', show);
f1b = get_distribution_fns(data.incd2020,   'lognorm', show);
f2  = get_distribution_fns(data.mort,       'lognorm', show);
f3  = get_distribution_fns(data.p_vulnpopn, 'beta',    show);
f4  = get_distribution_fns(data.p_vulnTB,   'beta',    show);
f5  = get_distribution_fns(data.propincd_ch,'beta', show);
f6  = get_distribution_fns(data.p_chpopn,  'beta',    show);
f7  = get_distribution_fns(data.p_adpopn,  'beta',    show);
f8  = get_distribution_fns(data.ch_notifs, 'lognorm', show);
f9 = get_distribution_fns(data.ART_covg,    'beta', show);
f10 = get_distribution_fns(data.HIV_prev,    'lognorm', show);


% lhd.fn = @(incd, mort, p_migrTB, p_migrpopn, p_LTBI) f1(incd) + f2(mort) + f3(p_migrTB) + f4(p_migrpopn) + f5(p_LTBI);
% lhd.fn = @(incd2010, incd2020, mort, p_migrTB, p_migrpopn, p_LTBI, nTPT2019) f1a(incd2010) + f1b(incd2020) + f2(mort) + f3(p_migrTB) + f4(p_migrpopn) + f5(p_LTBI) + f6(nTPT2019);
lhd.fn = @(incd2010, incd2020, mort, p_vulnpopn, p_vulnTB, p_chpopn, ch_notifs, ART_covg, HIV_prev) f1a(incd2010) + f1b(incd2020) + ...
                                                                                                    f2(mort) + f3(p_vulnpopn) + f4(p_vulnTB) ...
                                                                                                    + f6(p_chpopn) + f8(ch_notifs) ...
                                                                                                    + f9(ART_covg) + f10(HIV_prev);

save Model_setup;


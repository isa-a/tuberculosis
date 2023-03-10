clear all;

states   = {'U','Lf','Ls','I','Tx','Rlo','Rhi','R'};
gps.born = {'dom','for'};
gps.age  = {'ad','el'};

[i, s, d, lim] = get_addresses({states, gps.born, gps.age}, [], [], [], 0);
d = char(d);

% Include the auxiliaries
i.aux.inc  = i.nstates + [1:3];
i.aux.mort = i.nstates + 4;
i.nx = i.aux.mort(end);

s.infectious = [s.I];
s.prevalent  = [s.infectious, s.Tx];

% Selectors for the incidence
tmp = zeros(2,i.nstates); 
tmp(1,s.I) = 1;
tmp(2,intersect(s.I,s.for)) = 1;
tmp(3,intersect(s.I,s.el))  = 1;
agg.inc = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.I,:) = 1;
tmp(s.el, s.ad) = 0;
sel.inc = tmp - diag(diag(tmp));


% -- Natural history parameters -------------------------------------------
r.progression0  = 0.0826;
r.LTBI_stabil   = 0.872;
r.reactivation0 = 0.0006;

r.Tx            = 2;
r.default       = 0.01;

r.self_cure    = 1/6;
r.relapse      = [0.032 0.14 0.0015];
r.mu           = 1/83;                                                     % natural mortality
r.muTB         = 1/6;                                                      % TB related mortality
p.imm          = 0.8;                                                      % Reduced susceptibility conferred by previous infection


% -------------------------------------------------------------------------
% --- Name free parameters ------------------------------------------------

names = {'beta','gamma','p_birth','p_kLf','eld_relrisk','adu_mort'};
lgths =      [1,      1,       1,       1,            1,         1];

lim = 0; xi = [];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    xi.(names{ii}) = inds;
    lim = inds(end);
end

% Set their boundaries
bds = [];
bds(xi.beta,:)        = [0 30];
bds(xi.gamma,:)       = [0 50];
bds(xi.p_birth,:)     = [0 0.6];
bds(xi.p_kLf,:)       = [1 500];
bds(xi.eld_relrisk,:) = [1 500];
bds(xi.eld_relrisk,:) = [0 100];
bds(xi.adu_mort,:)    = [0 0.1];
prm.bounds = bds';

ref.i = i; ref.s = s; ref.xi = xi;
prm.p = p; prm.r = r; prm.agg = agg; prm.sel = sel;

% -------------------------------------------------------------------------
% --- Specify data --------------------------------------------------------

iso3 = 'AUS';

% Incidence and mortality
C = readtable('../Data/TB_burden_countries_2022-04-12.csv');
colnames = C.Properties.VariableNames;

cols = {'iso3','year','e_inc_100k_lo','e_inc_100k','e_inc_100k_hi','e_mort_100k_lo','e_mort_100k','e_mort_100k_hi'};
colnos = zeros(1,length(cols));
for ii = 1:length(colnos)
   colnos(ii) = find(strcmp(colnames,cols{ii})); 
end
C2  = C(:,colnos);
C2b = C2(C2.year==2019,:);

row = find(strcmp(C2b.iso3,iso3));
data.incd = C2b{row,3:5};
data.mort = C2b{row,6:8};

% --- Contribution of migrants
load ../Data/country_data.mat;
row = find(strcmp(tbl_migr.iso3,iso3));
if ~isempty(row)
    data.p_migrTB = tbl_migr{row,2}*[0.95 1 1.05];
else
    hm;
end

data.p_migrpopn = [0.27 0.30 0.33];
data.p_LTBI     = [0.15 0.2 0.25];

% --- Contribution of elderly
row = find(strcmp(tbl_old.ctrs2,iso3));
if ~isempty(row)
    data.p_eldTB = tbl_old{row,2}*[0.95 1 1.05];
else
    hm;
end
data.p_eldpopn = [0.21 0.23 0.25];


% --- Construct the calibration functions
show = 0;
f1 = get_distribution_fns(data.incd,       'lognorm', show);
f2 = get_distribution_fns(data.mort,       'lognorm', show);
f3 = get_distribution_fns(data.p_migrTB,   'beta', show);
f4 = get_distribution_fns(data.p_migrpopn, 'beta', show);
f5 = get_distribution_fns(data.p_LTBI,     'beta', show);
f6 = get_distribution_fns(data.p_eldTB,    'beta', show);
f7 = get_distribution_fns(data.p_eldpopn,  'beta', show);

lhd.fn = @(incd, mort, p_migrTB, p_migrpopn, p_LTBI, p_eldTB, p_eldpopn) f1(incd) + f2(mort) + f3(p_migrTB) + f4(p_migrpopn) + f5(p_LTBI) + f6(p_eldTB) + f7(p_eldpopn);

save Model_setup;


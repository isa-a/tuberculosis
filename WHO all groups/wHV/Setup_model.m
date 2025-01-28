clear all;

Get_references_Full;
Get_references_Reduced;                                                    % Without vaccination

% --- Define variables and ranges -----------------------------------------
names = {'r_beta','r_sym','r_cs','r_cs2','p_Dx','rf_mort_TB','r_ART_init','r_HIV_mort','rf_self_cure','p_HIV_relrate'};   % <--- Updated
lgths =        [2,      1,     1,      1,     1,           2,           1,           1,             1,              1];

xi = []; lim = 0;
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    xi.(names{ii}) = inds;
    lim = inds(end);
end
xi.calib = xi.r_HIV_mort;

bds = [];
bds(xi.r_beta,:)        = [0 30; 0 1];
bds(xi.r_sym,:)         = [0.1 10];
bds(xi.r_cs,:)          = [0.1 30];
bds(xi.r_cs2,:)         = [1 24];
bds(xi.p_Dx,:)          = [0 1];
bds(xi.rf_mort_TB,:)    = [0 2; 0 20];
bds(xi.r_ART_init,:)    = [0 10];
bds(xi.r_HIV_mort,:)    = [0 10];
bds(xi.rf_self_cure,:)  = [0.8 1.2];
bds(xi.p_HIV_relrate,:) = [1 100];
prm.bounds = bds';



% --- Define baseline parameter values ------------------------------------

% Natural history
r.progression0  = 0.0826*[1 10 10*0.4].*[1 1 0.5];
r.LTBI_stabil   =  0.872*[1 0 1];
r.reactivation0 = 0.0006*[1 100 100*0.4].*[1 1 0.5];
r.self_cure0    = 1/6*[1 0 1];
r.mort_TB0      = 1/6;
r.relapse       = [0.032 0.14 0.0015];
r.mort          = 1/66;
p.imm           = [0.8 0 0.8];

% Diagnosis stage
r.Dx           = 52;

% Treatment
load ../../data_AFR Tx_outcomes;
Txvec      = Tx_outcomes([1,2],:);
pcomplete  = Txvec(:,1)+Txvec(:,2);
r.Tx       = 2;

tmp        = r.Tx.*Txvec(:,3)./pcomplete;   
r.muTx     = tmp([1 2 2]);

tmp        = r.Tx.*Txvec(:,4)./pcomplete;
r.ltfu     = tmp([1 2 2]);

tmp        = Txvec(:,1)./pcomplete;                                        % Proportion success out of those completing
p.succ     = tmp([1 2 2]);

% Interventions
r.TPT       = 0;
p.TPT_PLHIV = 0;
r.ACF       = [0 0];
r.vacc    = 0;
r.waning  = 1/10;
p.VE      = [0 0.5 0];

% Bring them all together
prm.p = p; prm.r = r;
ref.i = i; ref.s = s; ref.d = d; ref.xi = xi;

ref0.i = i0; ref0.s = s0; ref0.d = d0; ref0.xi = xi;

% --- Get calibration targets ---------------------------------------------
load ../../data_AFR data;
load ../../data_Thembisa_AFR.mat;

data.HIV_prev = HIVprev_2019;
data.ART_covg = ARTcovg_2019;
data.psym     = [0.4 0.5 0.6];
prm.ART_start = ART_start;

% Extrapolate HIV incidence by 5 years
ys1 = HIV_incd(:,2);
xs1 = 1:length(ys1);
xs2 = 1:length(ys1)+30;
ys2 = max(interp1(xs1,ys1,xs2,'linear','extrap'),0);
prm.rHIV = ys2;

show = 0;
f1 = get_distribution_fns(data.inc2019_all, 'lognorm', show);
f2 = get_distribution_fns(data.inc2019_H1,  'lognorm', show);
f3 = get_distribution_fns(data.notifs(1)*[0.85 1 1.15], 'lognorm', show);
f4 = get_distribution_fns(data.ART_covg,    'beta', show);
f5 = get_distribution_fns(data.HIV_prev,    'lognorm', show);
f6 = get_distribution_fns(data.mort2019_H0, 'lognorm', show);
f7 = get_distribution_fns(data.mort2019_H1, 'lognorm', show);
f8 = get_distribution_fns(data.psym,     'beta', show);

lhd.fn  = @(inc_all, inc_h1, noti, ART_covg, HIV_prev, mort_H0, mort_H1, psym) f1(inc_all) + f2(inc_h1) + f3(noti) + f4(ART_covg) + f5(HIV_prev) + f6(mort_H0) + f7(mort_H1) + f8(psym);
lhd.sgn = -Inf;


lhd.fn1 = @(x) f1(x);
lhd.fn2 = @(x) f2(x);
lhd.fn3 = @(x) f3(x);
lhd.fn4 = @(x) f4(x);
lhd.fn5 = @(x) f5(x);
lhd.fn6 = @(x) f6(x);
lhd.fn7 = @(x) f7(x);

save Model_setup;


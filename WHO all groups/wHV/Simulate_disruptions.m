clear all; load model_fits2; load Model_setup.mat;

obj = @(x) get_objective(x, prm, ref0, sel0, agg0, gps0, lhd);

ix0 = size(xsto,1)/2; nx = 150; dx = round((size(xsto,1)-ix0)/nx);
xs = xsto(ix0:dx:end,:,1);


% --- Get the South Africa notifications
load ../../data_SSAfrica.mat;
% --- Get the Regional notifications
notif_rate = notifs;


% --- Lockdown dates
ldown   = 2020 + [4 4.5]/12;                                               % Start and end dates of lockdown
tend2   = 2031;                                                            % End date for simulation
betared = 0.15 + (0.85-0.15)*rand(size(xs,1),1);

% --- Disruption parameters
tmp = ones(1,length(2019:tend2));
% tmp([1:5]) = [1 0.75 0.8 1 1.1];
% fac = tmp;
vec0 = [1 0.7 0.9 0.9 1 0.94 1.01 1.05];
tmp1 = repmat(vec0,3,1);
tmp2 = linspace(vec0(end),1,12); tmp2(1) = [];
vec  = [tmp1(:)', tmp2];
fac  = ones(1,(tend2-2019)*12); fac(12+[1:length(vec)]) = vec;

% --- Do the simulations
opts = odeset('NonNegative',[1:i0.nstates],'Refine',64,'AbsTol',1e-10,'RelTol',1e-10);
opts2 = odeset('NonNegative',[1:i0.nstates],'Refine',64,'AbsTol',1e-14,'RelTol',1e-14); warning off;

r0 = r; p0 = p;
inct = [];
noti = [];

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    if mod(ii,mk) == 0; fprintf('%0.5g ',ii/mk); end
    
    [r,p] = alloc_parameters(xs(ii,:),r0,p0,xi);
    
    % Get the initial conditions
    [out, aux] = obj(xs(ii,:));
    init = aux.soln(end,1:end-1);
    
    
    % --- Simulate in absence of disruption -------------------------------
    p0 = p; r0 = r;
    r0.ART_init = 0;
    M0 = make_model(p0, r0, i0, s0, gps0);
    
    % Start of ART onwards
    pe = p; re = r;
    Me = make_model(pe, re, i0, s0, gps0);
    
    geq = @(t,in) goveqs_scaleup(t, in, M0, Me, [prm.ART_start 2019], i0, prm, sel0, agg0);
    
    [t0, soln0] = ode15s(geq, [2019:tend2], init, opts);
    allsol = soln0;
        
    
    % --- Now simulate disruption -----------------------------------------
    geq       = @(t,in) goveqs_scaleup_disruption(t, in, M0, Me, [prm.ART_start 2019], fac, 1, i0, prm, sel0, agg0);
    geq_ldown = @(t,in) goveqs_scaleup_disruption(t, in, M0, Me, [prm.ART_start 2019], fac, (1-betared(ii)), i0, prm, sel0, agg0);

    % Pre-disruption
    [ta, solna] = ode15s(geq, [2019 ldown(1)], init, opts);
    
    % During lockdown
    initb = solna(end,:);
    [tb, solnb] = ode15s(geq_ldown, [ldown(1) ldown(2)], initb, opts);
    
    % After lockdown
    initc = solnb(end,:);
    [tc, solnc] = ode15s(geq, [ldown(2) tend2], initc, opts);
    
    soln   = [solna; solnb(2:end,:); solnc(2:end,:)];
    t      = [ta; tb(2:end); tc(2:end)];
    soln1  = interp1(t, soln, t0);
    allsol = cat(3,allsol,soln1);
        

    % ---------------------------------------------------------------------
    % --- Get the outputs for incidence etc
    
    inct(:,ii,:) = sum(diff(allsol(:,i0.aux.inc,:),1),2);
    mort(:,ii,:) = sum(diff(allsol(:,i0.aux.mort,:),1),2);
    noti(:,ii,:) = sum(diff(allsol(:,i0.aux.noti,:),1),2);
       
end
fprintf('\n');

inc_pct  = permute(prctile(inct,[2.5,50,97.5],2)*1e5,[2,1,3]);             % Dims: 1.Lo/Md/Hi 2.Month 3.Scenario
mrt_pct  = permute(prctile(mort,[2.5,50,97.5],2)*1e5,[2,1,3]);
noti_pct = permute(prctile(noti,[2.5,50,97.5],2)*1e5,[2,1,3]);           

mat = squeeze(noti_pct(2,:,:));
figure; hold on;

plot(2019:(tend2-1), mat);
plot(2019:2022, notifs, 'g.-', 'MarkerSize', 24);

yl = ylim; yl(1) = 0; ylim(yl);
xlim([2019 2022]);

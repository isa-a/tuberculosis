clear all; load model_fits2; load Model_setup.mat;

obj = @(x) get_objective2(x, prm, ref, sel, agg, gps, lhd);

nx = 50; 
ix0 = 2e3; dx = round((size(xsto,1)-ix0)/nx);
xs = xsto(ix0:dx:end,:,1);


% --- Lockdown dates
ldown   = 2020 + [4 4.5]/12;                                               % Start and end dates of lockdown
tend1   = 2025;
tend2   = 2036;                                                            % End date for simulation
betared = 0.15 + (0.85-0.15)*rand(size(xs,1),1);

% --- Disruption parameters
vec0 = [1 0.7 0.9 0.9 1 0.94 1.01 1.05];                                   % Determined using simulate_disruption
tmp1 = repmat(vec0,3,1);
tmp2 = linspace(vec0(end),1,12); tmp2(1) = [];
vec  = [tmp1(:)', tmp2];
fac  = ones(1,(tend2-2019+5)*12); fac(12+[1:length(vec)]) = vec;

% --- Do the simulations
opts = odeset('NonNegative',[1:i.nstates],'Refine',64,'AbsTol',1e-10,'RelTol',1e-10);
opts2 = odeset('NonNegative',[1:i.nstates],'Refine',64,'AbsTol',1e-14,'RelTol',1e-14); warning off;

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
    ppre = p; rpre = r;
    rpre.ART_init = 0;
    Mpre = make_model(ppre, rpre, i, s, gps);
    
    % Start of ART onwards
    p0 = p; r0 = r;
    M0 = make_model(p0, r0, i, s, gps);
    
    geq = @(t,in) goveqs_scaleup(t, in, Mpre, M0, [prm.ART_start 2019], i, prm, sel, agg);
    
    [t0, soln0] = ode15s(geq, [2019:tend2], init, opts);
    allsol = soln0;
    

    % --- Now simulate pre-intervention -----------------------------------
    geq       = @(t,in) goveqs_scaleup_disruption(t, in, Mpre, M0, [prm.ART_start 2019], fac, 1, i, prm, sel, agg);
    geq_ldown = @(t,in) goveqs_scaleup_disruption(t, in, Mpre, M0, [prm.ART_start 2019], fac, (1-betared(ii)), i, prm, sel, agg);
    
    % Pre-disruption
    [ta, solna] = ode15s(geq, [2019 ldown(1)], init, opts);
    
    % During lockdown
    initb = solna(end,:);
    [tb, solnb] = ode15s(geq_ldown, [ldown(1) ldown(2)], initb, opts);
    
    % After lockdown
    initc = solnb(end,:);
    [tc, solnc] = ode15s(geq, [ldown(2) tend1], initc, opts);
    
    % Construct solution so far
    tmp   = [solna; solnb(2:end,:); solnc(2:end,:)];
    t     =    [ta;      tb(2:end);      tc(2:end)];
    soln  = interp1(t, tmp, 2019:tend1);


    initd = soln(end,:);
    % --- Modelling interventions -----------------------------------------
    
    % Improved diagnosis
    p1 = p; r1 = r;
    p1.Dx = 0.9;
    M1 = make_model(p1, r1, i, s, gps);

    % Accelerated case-finding, symptomatic only
    p2 = p1; r2 = r1;
    r2.ACF = [1, 0];
    M2 = make_model(p2, r2, i, s, gps);
    
    % Accelerated case-finding, asymptomatic as well
    p3 = p2; r3 = r2;
    r3.ACF = [1 1];
    M3 = make_model(p3, r3, i, s, gps);

    % Improve uptake of ART
    p4 = p3; r4 = r3;
    r4.ART_init = 2*r3.ART_init;
    M4 = make_model(p4, r4, i, s, gps);

    % Model them all 
    models = {M0, M1, M2, M3, M4};
    for mi = 1:length(models)
        geq = @(t,in) goveqs_scaleup(t, in, M0, models{mi}, tend1 + [0 3], i, prm, sel, agg);
        [td, solnd] = ode15s(geq, [tend1:tend2], initd, opts);
        allsol(:,:,mi+1) = [soln; solnd(2:end,:)];
    end
    

    inct(:,:,ii)     = squeeze(sum(diff(allsol(:,i.aux.inc,:),1),2));
    inct_hiv(:,:,ii) = squeeze(sum(diff(allsol(:,i.aux.inc([2,3]),:),1),2));
    mort(:,:,ii)     = squeeze(sum(diff(allsol(:,i.aux.mort,:),1),2));
    noti(:,:,ii)     = squeeze(sum(diff(allsol(:,i.aux.noti,:),1),2));
       
end
fprintf('\n');

inc_pct    = permute(prctile(inct,[2.5,50,97.5],3)*1e5,[3,1,2]);             % Dims: 1.Lo/Md/Hi 2.Month 3.Scenario
inchiv_pct = permute(prctile(inct_hiv,[2.5,50,97.5],3)*1e5,[3,1,2]);             % Dims: 1.Lo/Md/Hi 2.Month 3.Scenario
mrt_pct    = permute(prctile(mort,[2.5,50,97.5],3)*1e5,[3,1,2]);
noti_pct   = permute(prctile(noti,[2.5,50,97.5],3)*1e5,[3,1,2]);           

mat  = squeeze(inc_pct(2,:,:));
mat2 = squeeze(inchiv_pct(2,:,:));

figure; hold on;
plot(2019:(tend2-1),mat,'LineWidth',1.5);
plot(2019:(tend2-1),mat2,'--','LineWidth',1.5);
legend('No covid','Covid','Improved Dx','ACF (sympto)','ACF (all)','ART increase','location','SouthWest');
yl = ylim; yl(1) = 0; ylim(yl);



% save res_forward2;
% Plot_interventions;
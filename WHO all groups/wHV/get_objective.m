% v2: Version to include scale-up of TPT amongst PLHIV

function [out, aux, msg] = get_objective(x, prm, ref, sel, agg, gps, calfn)

r = prm.r; p = prm.p; i = ref.i; s = ref.s; xi = ref.xi;
mat = [prm.bounds(2,1:length(x))-x; x-prm.bounds(1,1:length(x))];

if min(mat(:)) < 0
    out = Inf*calfn.sgn;
    aux = [];
    msg = 1;
else
    
    [r,p] = alloc_parameters(x,r,p,xi);
    
    % --- Set up the necessary models -------------------------------------
    
    % Pre-ART conditions
    p0 = p; r0 = r;
    r0.ART_init = 0;
    M0 = make_model(p0, r0, i, s, gps);
      
    % Start of ART, but pre-TPT
    p1 = p; r1 = r;
    M1 = make_model(p1, r1, i, s, gps);
    
    % Start of ART, and with TPT in PLHIV
    p2 = p1; r2 = r1;
    p2.TPT_PLHIV = 0.58;
    M2 = make_model(p2, r2, i, s, gps);

    % --- Simulate the models
    
    % Equilibrium model
    init = zeros(1,i.nx); 
    seed = 1e-6;
    init(i.U.v0.h0) = (1-seed); 
    init(i.I.v0.h0) = seed;
    geq = @(t,in)goveqs_basis2(t, in, M0, i, prm, sel, agg);
    [t0, soln0] = ode15s(geq, [0:2e3], init, odeset('NonNegative',[1:i.nstates],'RelTol',1e-10,'AbsTol',1e-10));
    
    % HIV decline/ART scaleup model
    init = soln0(end,:);
    geq = @(t,in) goveqs_scaleup2D(t, in, M0, M1, M2, [prm.ART_start 2019; 2011 2019], i, prm, sel, agg);
    [t1, soln1] = ode15s(geq, [prm.ART_start:2019], init, odeset('NonNegative',[1:i.nstates],'RelTol',1e-10,'AbsTol',1e-10));
        
    sfin  = soln1(end,:);
    sdiff = diff(soln1,1);
    
    % --- Get the objectives ----------------------------------------------
    
    inc_all  = sum(sdiff(end,i.aux.inc))*1e5;
    inc_h1   = sum(sdiff(end,i.aux.inc([2,3])))*1e5;
    noti     = sum(sdiff(end,i.aux.noti))*1e5;
    ART_covg = sum(sfin(s.hart)/sum(sfin([s.h1,s.hart])));
    HIV_prev = sum(sfin([s.h1, s.hart]))/sum(sfin(1:i.nstates));
    mort_H0  = sdiff(end,i.aux.mort(1))*1e5;
    mort_H1  = sdiff(end,i.aux.mort(2))*1e5;
    psym     = sum(sfin(s.symptomatic))/sum(sfin(s.prevalent)); 
    
    if inc_all < 1
        out = Inf*calfn.sgn;
        aux = [];
        msg = 2;
    else
        out = calfn.fn(inc_all, inc_h1, noti, ART_covg, HIV_prev, mort_H0, mort_H1, psym);

        % --- Get additional outputs and package --------------------------
        aux.soln     = [soln1, t1];
        aux.allsol   = [soln0; soln1(2:end,:)];
        aux.inc_all  = inc_all;
        aux.inc_h1   = inc_h1;
        aux.noti     = noti;
        aux.ART_covg = ART_covg;
        aux.HIV_prev = HIV_prev;
        aux.mort_H0  = mort_H0;
        aux.mort_H1  = mort_H1;
%         aux.sym      = sym_sim;
        aux.M1       = M1;
        aux.sim      = [inc_all, inc_h1, noti, ART_covg, HIV_prev, mort_H0, mort_H1, psym];    
        msg = 0;
    end        
end

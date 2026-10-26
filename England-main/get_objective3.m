% v3: Now dropping vulnerable groups, and also tidying redundant terms
% relating to RR-TB

function [out, aux, msg] = get_objective2(x, ref, prm, gps, rin_vec, calfn)

i = ref.i; s = ref.s; xi = ref.xi;
p = prm.p; r = prm.r; sel = prm.sel; agg = prm.agg;

[p, r, prm] = allocate_parameters(x, p, r, xi, prm);

% keyboard;

tmp1  = [prm.bounds; x];
tmp2  = diff(tmp1([1,3,2],:),[],1);
cond1 = min(tmp2(:))<0;

if cond1
    out = -Inf;
    aux = NaN;
    msg = 0;
else
    init = zeros(1,i.nx);
    seed = 1e-5;
    init(i.U.dom)    = 1 - seed;
    init(i.Irec.dom) = seed;
    
    % Equlibrium model
    p0 = p; r0 = r; 
    p0.LTBIdec = 0;
    p0.prev_in_migr = 0.003;
    r0.gamma   = r.gamma_2015;
    M0 = make_model(p0, r0, i, s, gps, prm.contmat);
    
    % >2015: scaleup of TPT 
    p1 = p; r1 = r; 
    r1.TPT = [0 r.TPT2020rec 0];
    r1.gamma = r.gamma_2015;
    M1 = make_model(p1, r1, i, s, gps, prm.contmat);
    
    % >2010: increase in case-finding, and decline of beta
    p2 = p; r2 = r; 
    r2.gamma = r.gamma_2020;
    p2.prev_in_migr = 0;
    % r2.gamma = r.gamma_2015;
    M2 = make_model(p2, r2, i, s, gps, prm.contmat);

    % --- Now simulate them all

    geq0 = @(t,in) goveqs_basis3(t, in, i, s, M0, rin_vec, agg, sel, r0, p0, true);
    [t0, soln0] = ode15s(geq0, [0:5e3], init, odeset('NonNegative', 1:i.nstates));

    % Increased TPT and case-finding
    geq1 = @(t,in) goveqs_scaleup2D(t, in, M0, M1, M2, rin_vec, [2015 2020; 2010 2011], i, s, p2, r2, prm, sel, agg, false);
    [t1, soln1] = ode15s(geq1, [2010:2021], soln0(end,:), odeset('NonNegative', 1:i.nstates));
    
    if sum(soln1(end,1:i.nstates),2)>5
        keyboard;
    end

    dsol = diff(soln1,[],1);
    sfin = soln1(end,:);
    pops = sum(soln1(:,1:i.nstates),2);

    incd     = dsol(:,i.aux.inc)./pops(1:end-1)*1e5;
    incd2010 = incd(t1==2010,1);
    incd2020 = incd(t1==2020,1);
    p_migrTB = incd(end,2)/incd(end,1);

    mort     = dsol(:,i.aux.mort)./pops(1:end-1)*1e5;
    mort2020 = mort(t1==2020);
    
    p_migrpopn       = sum(sfin(s.migr))/pops(end);                        % Proportion of popn being migrants
    p_incd_recentinf = incd(end,3)/incd(end,1);                            % Proportion of incidence coming from recent infection 

    p_LTBI_inmigr = sum(sfin(intersect(s.migr_rect,[s.Lf, s.Ls])))/sum(sfin(s.migr_rect));

    % Number initiating TPT in 2019
    n_TPT2019  = dsol(end,i.aux.nTPT)*1e5/pops(end);
    

    % --- Also, do cohort simulation
    r3 = r2; p3 = p2;
    %r3.beta = 0;
    r3in_vec = 0*rin_vec;

    % Set up a cohort of recently-arrive migrants
    M3 = make_model(p3, r3, i, s, gps, prm.contmat);
    init = zeros(1, i.nx);
    init(1:i.nstates) = M3.migrentries';

    geq0 = @(t,in) goveqs_basis3(t, in, i, s, M3, r3in_vec, agg, sel, r3, p3, true);
    [coh_t, coh_soln] = ode15s(geq0, [0:10], init, odeset('NonNegative', 1:i.nstates));
    
    migr_incd_5yr = coh_soln((coh_t==5),i.aux.inc(1))*100;


    % --- Evaluate the posterior density 
    if incd(end,1) < 0.1
        out = -Inf;
        aux = NaN;
        msg = 1;
    else
        sim = [incd2010, incd2020, mort2020, p_migrTB, p_migrpopn, p_incd_recentinf, migr_incd_5yr];
        out = calfn.fn(sim);
        msg = 2;

        aux.sim        = sim;
        aux.init       = soln1(end,:);
        aux.soln       = soln1;
        aux.incdt      = incd;
        aux.nTPT       = n_TPT2019;
    end
end
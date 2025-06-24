function [out, aux, msg] = get_objective2(x, ref, prm, gps, contmat, calfn)

i = ref.i; s = ref.s; xi = ref.xi;
p = prm.p; r = prm.r; sel = prm.sel; agg = prm.agg;

[p,r,prm] = allocate_parameters(x, p, r, prm, xi);

% keyboard;

tmp1  = [prm.bounds; x];
tmp2  = diff(tmp1([1,3,2],:),[],1);
cond1 = min(tmp2(:))<0;

if cond1
    out = -Inf;
    aux = NaN;
    msg = 0;
else

    %Set up models

    % Equlibrium model, pre ART, pre HIV
    p0 = p; r0 = r; 
    p0.betadec = 0;
    r0.gamma   = r.gamma_2015;
    p0.relbeta = 0; r0.RR_acqu = 0;
    r0.ART_init = 0;
    M0 = make_model(p0, r0, i, s, gps, prm.contmat);

    % HIV starts but no ART 
    p1 = p; r1 = r; 
    r1.TPT = [0 r.TPT2020rec 0 0];
    r1.gamma = r.gamma_2015;
    r0.ART_init = 0;
    M1 = make_model(p1, r1, i, s, gps, prm.contmat);
     
    % >HIV, and ART scaleup
    p2 = p; r2 = r; 
    r2.gamma = r.gamma_2020;
    M2 = make_model(p2, r2, i, s, gps, prm.contmat);
    
    % --- Now simulate them all

    % Eq model
    init = zeros(1,i.nx);
    seed = 1e-5;
    init(i.U.ad.dom.neg)       = 1 - seed;
%     init(i.I.ad.dom.neg)    = seed;
    init(i.Irec.ad.dom.neg) = seed;
    options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10, 'NonNegative', 1:i.nstates);
    geq0 = @(t,in) goveqs_basisnonHIV(t, in, i, s, M0, agg, sel, r0, p0);
    [t0, soln0] = ode15s(geq0, [0:5e3], init, options);

    % HIV decline/ART scaleup model
    init = soln0(end, :);
    geq1 = @(t,in) goveqs_scaleup2D(t, in, M0, M1, M2, [prm.ART_start 2020; prm.ART_start 2020], i, s, p2, r2, prm, sel, agg);
    [t1, soln1] = ode15s(geq1, [prm.ART_start:2020], init, options);
    
    dsol   = diff(soln1,[],1);
    sfin   = soln1(end,:);
    
    %-- Record --------------------------------------------------------
    incd2010   = dsol(t1==2010,i.aux.inc(1))*1e5;
    incd2020   = dsol(end,i.aux.inc(1))*1e5;
    incd       = dsol(end,i.aux.inc)*1e5;
    mort       = dsol(end,i.aux.mort)*1e5;

    % Proportion of population being vulnerable
    %p_vulnpopn = sum(sfin(s.vuln))/sum(sfin(1:i.nstates));

    % Contribution of vulnerable people to overall incidence
    %p_vulnTB   = incd(2)/incd(1);

    % Number initiating TPT in 2019
    n_TPT2019  = dsol(end,i.aux.nTPT)*1e5;

    % Contribution of kids to overall incidence
%     propincd_ch =  incd(2)/incd(1);

    % Proportion of population thats kids
%     p_chpopn = sum(sfin(s.ch))/sum(sfin(1:i.nstates));

    % Proportion of population thats adults
%     p_adpopn = sum(sfin(s.ad))/sum(sfin(1:i.nstates));

    % Notifications
%     ch_notifs = dsol(end,i.aux.ch_notifs)*1e5;

    % ART coverage
    ART_covg = sum(sfin(s.art)/sum(sfin([s.pos,s.art])));

    % HIV prevalence
    HIV_prev = sum(sfin([s.pos, s.art]))/sum(sfin(1:i.nstates));

    p_incd_recentinf = incd(end,5)/incd(end,1);

    if incd(1) > 0.1
        % keyboard;
        out  = calfn.fn(incd2010, incd2020, mort, ART_covg, HIV_prev, p_incd_recentinf);
        aux.soln       = soln1;
        msg            = 2;
        aux.incd       = dsol(find(t1==2010):end,i.aux.inc(1))*1e5;
        aux.incd2010   = incd2010;
        aux.incd2020   = incd2020;
        aux.mort       = mort;
        aux.nTPT       = n_TPT2019;
%         aux.propincd_ch = propincd_ch;
%         aux.chpopn     = p_chpopn;
%         aux.adpopn     = p_adpopn;
%         aux.ch_notifs  = ch_notifs;
        aux.ART_covg   = ART_covg;
        aux.HIV_prev   = HIV_prev;
        aux.p_incd_recentinf = p_incd_recentinf;
        aux.sim        = [incd2010, incd2020, mort, ART_covg, HIV_prev, p_incd_recentinf];
    else
        % keyboard;
        out = -Inf;
        aux = NaN;
        msg = 1;
        % disp('incd <= 0.1');
        % disp(incd);
        % aux.soln0      = soln0;
        % aux.soln0b     = soln0b;
        % aux.soln       = soln1;
        % aux.incd       = dsol(find(t1==2010):end,i.aux.inc(1))*1e5;
        % aux.incd2010   = incd2010;
        % aux.incd2020   = incd2020;
        % aux.incdRR2020 = incdRR2020;
        % aux.mort       = mort;
        % aux.p_migrTB   = p_migrTB;
        % aux.p_migrpopn = p_migrpopn;
        % aux.p_LTBI     = p_LTBI;
        % aux.p_migrect  = sum(sfin(s.migr_rect))/sum(sfin(1:i.nstates));
        % aux.nTPT       = n_TPT2019;
        % aux.M0         = M0;
    end
    
    %keyboard;
end
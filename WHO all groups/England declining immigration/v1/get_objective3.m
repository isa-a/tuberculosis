% v3: Now dropping vulnerable groups, and also tidying redundant terms
% relating to RR-TB

function [out, aux, msg] = get_objective2(x, ref, prm, gps, contmat, rin_vec, calfn)

i = ref.i; s = ref.s; xi = ref.xi;
p = prm.p; r = prm.r; sel = prm.sel; agg = prm.agg;

[p,r] = allocate_parameters(x, p, r, xi, prm.scaling);

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
    init(i.U.ad.dom)       = 1 - seed;
    % init(i.U.ad.migr_rect) = p.migrect_popn;
    init(i.Irec.ad.dom.ds) = seed;
    
    % Equlibrium model, without RR-TB
    p0 = p; r0 = r; 
    p0.betadec = 0;
    r0.gamma   = r.gamma_2015;
    p0.relbeta = 0; r0.RR_acqu = 0;
    M0 = make_model(p0, r0, i, s, gps, contmat);

    % Introduction of RR-TB from 1970
%     p0b = p; r0b = r; 
%     p0b.betadec = 0;
%     r0b.gamma   = r.gamma_2015;
%     M0b = make_model(p0b, r0b, i, s, gps, contmat);
    
    % >2015: scaleup of TPT 
    p1 = p; r1 = r; 
    r1.TPT = [0 r.TPT2020rec 0];
    r1.gamma = r.gamma_2015;
    M1 = make_modelnonEQ(p1, r1, i, s, gps, contmat);
    
    % >2010: increase in case-finding
    p2 = p; r2 = r; 
    % r2.gamma = r.gamma_2020;
    r2.gamma = r.gamma_2015;
    M2 = make_modelnonEQ(p2, r2, i, s, gps, contmat);

    % --- Now simulate them all

    geq0 = @(t,in) goveqs_basis3(t, in, i, s, M0, rin_vec, agg, sel, r0, p0, true);
    [t0, soln0] = ode15s(geq0, [0:5e3], init, odeset('NonNegative', 1:i.nstates));

    % Emergence of RR-TB
    %geq0b = @(t,in) goveqs_basis3(t, in, i, s, M0b, agg, sel, r0b, p0b);
    %[t0b, soln0b] = ode15s(geq0b, [1970:2010], soln0(end,:), odeset('NonNegative',1:i.nstates));

    % Increased TPT and case-finding
    geq1 = @(t,in) goveqs_scaleup2D(t, in, M0, M1, M2, rin_vec, [2015 2020; 2010 2020], i, s, p2, r2, prm, sel, agg, false);
    [t1, soln1] = ode15s(geq1, [2010:2030], soln0(end,:), odeset('NonNegative', 1:i.nstates));
    
%     allsol = [soln0; soln1(2:end,:)];
%     allt   = [t0; t1(2:end)];
    dsol   = diff(soln1,[],1);
    sfin   = soln1(end,:);
    popfin = sum(sfin(1:i.nstates));

    incd2010   = dsol(t1==2010,i.aux.inc(1))*1e5/popfin;
    incd2020   = dsol(end,i.aux.inc(1))*1e5/popfin;
    %incdRR2020 = dsol(end,i.aux.inc(3))*1e5;
    incd       = dsol(end,i.aux.inc)*1e5/popfin;
    mort       = dsol(end,i.aux.mort)*1e5/popfin;
    p_migrTB   = incd(2)/incd(1);
    
    %p_LTBI     = sum(sfin(intersect(s.migr_rect,[s.Lf, s.Ls])))/sum(sfin(s.migr_rect));
    p_LTBI_inmigr = sum(sfin(intersect(s.migr_rect,[s.Lf, s.Ls])))/sum(sfin(s.migr_rect));
    p_migrpopn    = sum(sfin(s.migr))/popfin;

    % % Proportion of population being vulnerable
    % p_vulnpopn = sum(sfin(s.vuln))/sum(sfin(1:i.nstates));
    % 
    % % Contribution of vulnerable people to overall incidence
    % p_vulnTB   = incd(3)/incd(1);

    % Number initiating TPT in 2019
    n_TPT2019  = dsol(end,i.aux.nTPT)*1e5/popfin;

    % Proportion of population thats kids
    p_chpopn = sum(sfin(s.ch))/popfin;

    % Incidence in children; 2020
    incd_ch2020 = dsol(end,i.aux.inc(3))/sum(sfin(s.ch))*1e5;

    % % Proportion of population thats adults
    % p_adpopn = sum(sfin(s.ad))/sum(sfin(1:i.nstates));

    % Notifications
    ch_notifs = dsol(end,i.aux.ch_notifs)/popfin*1e5;

    % Proportion of incidence coming from recent infection
    p_incd_recentinf = incd(4)/incd(1);
    
    if incd(1) > 0.1
        % out  = calfn.fn(incd2010, incd2020, mort, p_migrTB, p_LTBI_inmigr, p_vulnpopn, p_vulnTB, incd_ch2020, p_chpopn, p_adpopn, ch_notifs);
        out  = calfn.fn(incd2010, incd2020, mort, p_migrTB, p_migrpopn, p_LTBI_inmigr, p_chpopn, ch_notifs, p_incd_recentinf);
        aux.soln0      = soln0;
        aux.soln       = soln1;
        msg            = 2;
        aux.incd       = dsol(find(t1==2010):end,i.aux.inc(1))/popfin*1e5;
        aux.incd2010   = incd2010;
        aux.incd2020   = incd2020;
        aux.mort       = mort;
        aux.p_migrTB   = p_migrTB;
        aux.p_migrpopn = p_migrpopn;
        aux.p_LTBI_inmigr = p_LTBI_inmigr;
        % aux.p_vulnpopn = p_vulnpopn;
        % aux.p_vulnTB   = p_vulnTB;
        aux.p_migrect  = sum(sfin(s.migr_rect))/sum(sfin(1:i.nstates));
        aux.nTPT       = n_TPT2019;
        aux.incd_ch    = incd_ch2020;
        aux.chpopn     = p_chpopn;
        % aux.adpopn     = p_adpopn;
        aux.ch_notifs  = ch_notifs;
        aux.p_incd_recentinf = p_incd_recentinf;
        aux.sim        = [incd2010, incd2020, mort, p_migrTB, p_migrpopn, p_LTBI_inmigr, incd_ch2020, p_chpopn, ch_notifs, p_incd_recentinf];
    else
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
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
    init = zeros(1,i.nx);
    seed = 1e-5;
    init(i.U.ad.dom.neg)       = 1 - seed;
    %init(i.U.ad.migr_rect) = p.migrect_popn;
    init(i.I.ad.dom.ds.neg)    = seed;
    
    % Equlibrium model, without RR-TB
    p0 = p; r0 = r; 
    p0.betadec = 0;
    r0.gamma   = r.gamma_2015;
    p0.relbeta = 0; r0.RR_acqu = 0;
    M0 = make_model(p0, r0, i, s, gps, prm.contmat);

    %Introduction of HIV from 1990 to 2010
    p1 = p; r1 = r; 
    r1.HIV = r.HIVincdpeak;
    M1 = make_model(p1, r1, i, s, gps, prm.contmat);

    % Introduction of RR-TB from 1970
%     p0b = p; r0b = r; 
%     p0b.betadec = 0;
%     r0b.gamma   = r.gamma_2015;
%     M0b = make_model(p0b, r0b, i, s, gps, contmat);
    
    % >2015: scaleup of TPT 
    p2 = p; r2 = r; 
    r2.TPT = [0 r.TPT2020rec 0 0];
    r2.gamma = r.gamma_2015;
    M2 = make_model(p2, r2, i, s, gps, prm.contmat);
    
    % >2010: increase in case-finding
    p3 = p; r3 = r; 
    r3.gamma = r.gamma_2020;
    M3 = make_model(p3, r3, i, s, gps, prm.contmat);

    % >2010: scaleup of ART 
    p4 = p; r4 = r; 
    %r3.TPT = [0 r.TPT2020rec 0 0];
    r4.ART = r.ARTnow;
    M4 = make_model(p4, r4, i, s, gps, prm.contmat);

    %decay of HIV peak from 2010 until 2020
    p5 = p; r5 = r; 
    r5.HIV = r.HIVincdnow;
    M5 = make_model(p5, r5, i, s, gps, prm.contmat);


    % --- Now simulate them all
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'NonNegative', 1:i.nstates);


    geq0 = @(t,in) goveqs_basis3(t, in, i, s, M0, agg, sel, r0, p0);
    [t0, soln0] = ode15s(geq0, [0:5e3], init, options);

    % Emergence of RR-TB
    %geq0b = @(t,in) goveqs_basis3(t, in, i, s, M0b, agg, sel, r0b, p0b);
    %[t0b, soln0b] = ode15s(geq0b, [1970:2010], soln0(end,:), odeset('NonNegative',1:i.nstates));

     % HIV scaled up to peak year
    geq0a = @(t,in) goveqs_scaleup1D(t, in, M0, M1, [1990:2010], i, s, p1, r1, prm, sel, agg);
    [t0a, soln0a] = ode15s(geq0a, [1990:2010], soln0(end,:), odeset('NonNegative',1:i.nstates));

    % Increased TPT and case-finding and  HIV decay
    geq1 = @(t,in) goveqs_scaleup2D(t, in, M1, M2, M3, M4, M5, [2015 2023; 2010 2023; 2010 2023; 2010 2023], i, s, p5, r5, prm, sel, agg);
    [t1, soln1] = ode15s(geq1, [2010:2023], soln0a(end,:), options);
    
%     allsol = [soln0; soln1(2:end,:)];
%     allt   = [t0; t1(2:end)];
    dsol   = diff(soln1,[],1);
    sfin   = soln1(end,:);
    
    incd2010   = dsol(t1==2010,i.aux.inc(1))*1e5;
    incd2020   = dsol(end,i.aux.inc(1))*1e5;
    %incdRR2020 = dsol(end,i.aux.inc(3))*1e5;
    incd       = dsol(end,i.aux.inc)*1e5;
    mort       = dsol(end,i.aux.mort)*1e5;
    %p_migrTB   = incd(2)/incd(1);
    
    %p_LTBI     = sum(sfin(intersect(s.migr_rect,[s.Lf, s.Ls])))/sum(sfin(s.migr_rect));
%     p_LTBI_inmigr = p.LTBI_in_migrch*p.ch_in_migr + p.LTBI_in_migrad*(1-p.ch_in_migr);
%     p_migrpopn    = sum(sfin(s.migr))/sum(sfin(1:i.nstates));

    % Proportion of population being vulnerable
    p_vulnpopn = sum(sfin(s.vuln))/sum(sfin(1:i.nstates));

    % Contribution of vulnerable people to overall incidence
    p_vulnTB   = incd(3)/incd(1);

    % Number initiating TPT in 2019
    n_TPT2019  = dsol(end,i.aux.nTPT)*1e5;

        % Incidence in children; 2020
%     incd_ch2020 = dsol(end,i.aux.inc(4))*1e5/p_chpopn;
%   propincd_ch = dsol(end,i.aux.inc(4)) / dsol(end,i.aux.inc(1));

    % Contribution of kids to overall incidence
    propincd_ch =  incd(4)/incd(1);

    % Proportion of population thats kids
    p_chpopn = sum(sfin(s.ch))/sum(sfin(1:i.nstates));

    % Proportion of population thats adults
    p_adpopn = sum(sfin(s.ad))/sum(sfin(1:i.nstates));

    % Notifications
    ch_notifs = dsol(end,i.aux.ch_notifs)*1e5;

    HIVincdpeak = dsol(t1==2010,i.aux.HIVinc)*1e5;
    HIVincdnow  = dsol(t1==2023,i.aux.HIVinc)*1e5;

    % incidence of tb amongst HIV positive and on ART
    incdTBHIV = dsol(end,i.aux.inc(5))*1e5;

    %coverage of ART
    ARTcov = sum(sfin(s.art))/sum(sfin([s.pos, s.art]));


    if incd > 0.1
        out  = calfn.fn(incd2010, incd2020, mort, p_vulnpopn, p_vulnTB, p_chpopn, ch_notifs, HIVincdpeak, HIVincdnow, incdTBHIV);
        aux.soln       = soln1;
        msg            = 2;
        aux.incd       = dsol(find(t1==2010):end,i.aux.inc(1))*1e5;
        aux.incd2010   = incd2010;
        aux.incd2020   = incd2020;
        aux.mort       = mort;
%         aux.p_migrTB   = p_migrTB;
%         aux.p_migrpopn = p_migrpopn;
%         aux.p_LTBI_inmigr = p_LTBI_inmigr;
        aux.p_vulnpopn = p_vulnpopn;
        aux.p_vulnTB   = p_vulnTB;
%         aux.p_migrect  = sum(sfin(s.migr_rect))/sum(sfin(1:i.nstates));
%         aux.p_migrect_ch  = sum(sfin(intersect(s.migr_rect,[s.ch])))/sum(sfin(s.migr_rect));
        aux.nTPT       = n_TPT2019;
        aux.propincd_ch = propincd_ch;
        aux.chpopn     = p_chpopn;
        aux.adpopn     = p_adpopn;
        aux.ch_notifs  = ch_notifs;
%         aux.vuln_prev  = vuln_prev;
%         aux.vuln_relrisk    = vuln_relrisk;
        aux.HIVincdpeak = HIVincdpeak;
        aux.HIVincdnow  = HIVincdnow;
        aux.incdTBHIV = incdTBHIV;
        aux.ARTcov    = ARTcov;
        aux.sim        = [incd2010, incd2020, mort, p_vulnpopn, p_vulnTB, propincd_ch, p_chpopn, p_adpopn, ch_notifs, HIVincdpeak, HIVincdnow, incdTBHIV, ARTcov];
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
function M = make_model(p, r, i, s, gps)

% --- Get the linear rates ------------------------------------------------
m  = zeros(i.nstates);
m2 = zeros(i.nstates);
m3 = zeros(i.nstates);

for iv = 1:length(gps.vacc)
    vacc = gps.vacc{iv};

    for ih = 1:length(gps.hiv)
        hiv = gps.hiv{ih};
        getst = @(st) i.(st).(vacc).(hiv);

        Lf  = getst('Lf');
        Ls  = getst('Ls');
        Lfp = getst('Lfp');
        Lsp = getst('Lsp');
        I   = getst('I');
        Isc = getst('Isc');
        E   = getst('E');
        Dx  = getst('Dx');
        Tx  = getst('Tx');
        Rlo = getst('Rlo');
        Rhi = getst('Rhi');
        R   = getst('R');

        % --- Fast progression and LTBI stabilisation
        source  = Lf;
        destins =                                   [Isc,                Ls];
        rates   = [r.progression(ih)*(1-p.VE(2)*(iv==2)), r.LTBI_stabil(ih)];
        m(destins, source) = m(destins, source) + rates';

        % --- Reactivation
        source = Ls;
        destin = Isc;
        rate   = r.reactivation(ih)*(1-p.VE(2)*(iv==2));
        m(destin, source) = m(destin, source) + rate;

        % --- Fast progression and LTBI stabilisation those who are under preventive therapy
        source  = Lfp;
        destins =                                           [Isc,               Lsp];
        rates   = [r.progression(ih)*(1-0.6)*(1-p.VE(2)*(iv==2)), r.LTBI_stabil(ih)];

        %keyboard;
        m(destins, source) = m(destins, source) + rates';

        % --- Reactivation from those who are under preventive therapy
        source = Lsp;
        destin = Isc;
        rate   = r.reactivation(ih)*(1-0.6)*(1-p.VE(2)*(iv==2));
        m(destin, source) = m(destin, source) + rate;

        % --- Sub-clinical to clinical disease
        source = Isc;
        destin = I;
        rate   = r.sym;
        m(destin, source) = m(destin, source) + rate;

        % --- Primary careseeking
        source = I;
        destin = Dx;
        rate   = r.cs;
        m2(destin, source) = m2(destin, source) + rate;

        % --- Secondary careseeking
        source = E;
        destin = Dx;
        rate   = r.cs2;
        m2(destin, source) = m2(destin, source) + rate;

        % --- Diagnosis outcomes
        source  = Dx;
        destins =        [Tx,      E];
        rates   = r.Dx*[p.Dx, 1-p.Dx];
        m(destins, source) = m(destins, source) + rates';

        % --- Treatment outcomes
        source  = Tx;
        destins =             [Rlo,                    E,         Rhi];
        rates   = [r.Tx*p.succ(ih),  r.Tx*(1-p.succ(ih)),  r.ltfu(ih)];
        m(destins, source) = m(destins, source) + rates';

        % --- Relapse
        sources = [Rlo, Rhi, R];
        destin  = Isc;
        rates   = r.relapse;
        m(destin, sources) = m(destin, sources) + rates;

        sources = [Rlo, Rhi];
        destin  = R;
        rates   = 0.5;
        m(destin, sources) = m(destin, sources) + rates;

        % --- Self cure
        sources = intersect(intersect(s.infectious,s.(hiv)),s.(vacc));
        destin  = Rhi;
        rates   = r.self_cure(ih);
        m(destin, sources) = m(destin, sources) + rates;

        % ---------------------------------------------------------------------
        % --- Interventions

        % Preventive therapy
        source = Ls; destin = Lsp; rate = r.TPT;
        m(destin, source) = m(destin, source) + rate;

        source = Lf; destin = Lfp; rate = r.TPT;
        m(destin, source) = m(destin, source) + rate;

        % ACF
        sources = [I, E, Dx];
        destin  = Tx;
        rate    = r.ACF(1);
        m(destin, sources) = m(destin, sources) + rate;

        sources = Isc;
        destin  = Tx;
        rate    = r.ACF(2);
        m(destin, sources) = m(destin, sources) + rate;

    end
end

% --- HIV acquisition
sources = s.h0;
destins = s.h1;
inds    = sub2ind(size(m), destins, sources);
rates   = 1;
m3(inds) = m3(inds) + rates;

% --- ART initiation
sources = s.h1;
destins = s.hart;
inds    = sub2ind(size(m), destins, sources);
rates   = r.ART_init;
m(inds) = m(inds) + rates;

% People being initiated on TPT as same time as ART
sources = intersect(s.h1,  [s.Lf,  s.Ls]);
destins = intersect(s.hart,[s.Lfp, s.Lsp]);
inds    = sub2ind(size(m), destins, sources);
rates   = r.ART_init*p.TPT_PLHIV;
m(inds) = m(inds) + rates;

% People missing TPT while initiating ART
sources = intersect(s.h1,  [s.Lf, s.Ls]);
destins = intersect(s.hart,[s.Lf, s.Ls]);
inds    = sub2ind(size(m), destins, sources);
rates   = r.ART_init*(1-p.TPT_PLHIV);
m(inds) = m(inds) + rates;

% Everyone else
sources = setdiff(s.h1,  [s.Lf, s.Ls]);
destins = setdiff(s.hart,[s.Lf, s.Ls]);
inds    = sub2ind(size(m), destins, sources);
rates   = r.ART_init;
m(inds) = m(inds) + rates;

if length(gps.vacc)>1
    % --- Uptake of vaccination
    sources = s.v0;
    destins = s.v1;
    rate    = r.vacc;
    inds    = sub2ind(size(m), destins, sources);
    m(inds) = m(inds) + rate;

    % --- Waning of vaccination
    sources = s.v1;
    destins = s.vw;
    rate    = r.waning;
    inds    = sub2ind(size(m), destins, sources);
    m(inds) = m(inds) + rate;
end

% --- Bring them together
M.lin    = sparse(m - diag(sum(m,1)));
M.Dxlin  = sparse(m2 - diag(sum(m2,1)));
M.linHIV = sparse(m3 - diag(sum(m3,1)));


% --- Get the nonlinear rates ---------------------------------------------

% --- Allocating transitions
m = zeros(i.nstates);
for iv = 1:length(gps.vacc)
    vacc = gps.vacc{iv};
    for ih = 1:length(gps.hiv)
        hiv = gps.hiv{ih};

        susinds = intersect(intersect([s.U, s.Lf, s.Ls, s.Lfp, s.Lsp, s.Rlo, s.Rhi, s.R],s.(hiv)),s.(vacc));
        imminds = setdiff(susinds, s.U);

        sources = susinds;
        destin  = i.Lf.(vacc).(hiv);
        m(destin, sources) = (1 - p.VE(1)*(iv==2));

        sources = intersect(intersect([s.Lfp, s.Lsp],s.(hiv)),s.(vacc));
        destin  = i.Lfp.(vacc).(hiv);
        m(destin, sources) = (1 - p.VE(1)*(iv==2));

        % Adjust for any immune protection
        m(:,imminds) = m(:,imminds)*(1-p.imm(ih));
    end
end
M.nlin = sparse(m - diag(sum(m,1)));


% --- Getting force-of-infection
m = zeros(1,i.nstates);
m(intersect(s.infectious,s.h0))          = r.beta(1);
m(intersect(s.infectious,[s.h1,s.hart])) = r.beta(2);
m(:,s.Isc) = m(:,s.Isc)*0.75;                                              % Subclinical TB have 75%transmission potential as clinical TB
M.lambda = sparse(m);


% --- Get the mortality rates
m = zeros(i.nstates,3);
m(:,1)            = r.mort;
m(s.h1,1)         = r.HIV_mort;                                             % HIV+ve, no TB
inds = intersect(s.h0,s.symptomatic);
m(inds,2) = r.mort_TB(1);                                                   % HIV-ve, untreated TB
inds = intersect([s.h1,s.hart],s.symptomatic);
m(inds,3) = r.mort_TB(2);                                                   % HIV+ve, untreated TB

inds = intersect(s.h0,s.Tx);
m(inds,2) = r.muTx(1);                                                      % HIV-ve, TB on treatment

inds = intersect([s.h1,s.hart],s.Tx);
m(inds,2) = r.muTx(2);                                                      % HIV+ve, TB on treatment


M.mortvec = m;
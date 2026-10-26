function M = make_model(p,r,i,s,gps,contmat)

m     = zeros(i.nstates);
mLTBI = zeros(i.nstates);

for ib = 1:length(gps.born)
    born = gps.born{ib};
    geti = @(st) i.(st).(born);

    Lf_imp = geti('Lf_imp');
    Lf     = geti('Lf');
    Ls     = geti('Ls');
    Pf_imp = geti('Pf_imp');
    Pf     = geti('Pf');
    Ps     = geti('Ps');
    Irem   = geti('Irem');
    Irec   = geti('Irec');
    Tx     = geti('Tx');
    Rlo    = geti('Rlo');
    Rhi    = geti('Rhi');
    R      = geti('R');

    % --- Outcomes of 'fast' latent
    source  = Lf;
    destins =            [Ls,              Irec];
    rates   = [r.LTBI_stabil, r.progression(ib)];
    m(destins, source) = m(destins, source) + rates';

    source  = Lf_imp;
    destins =            [Ls,              Irem];
    rates   = [r.LTBI_stabil, r.progression(ib)];
    m(destins, source) = m(destins, source) + rates';

    source  = Pf;
    destins =            [Ps,                           Irec];
    rates   = [r.LTBI_stabil, r.progression(ib)*(1-p.TPTeff)];
    m(destins, source) = m(destins, source) + rates';

    source  = Pf_imp;
    destins =            [Ps,                           Irem];
    rates   = [r.LTBI_stabil, r.progression(ib)*(1-p.TPTeff)];
    m(destins, source) = m(destins, source) + rates';

    % --- Outcomes of 'slow' latent
    source  = Ls;
    destin  = Irem;
    rate    = r.reactivation(ib);
    mLTBI(destin, source) = mLTBI(destin, source) + rate;

    source  = Ps;
    destin  = Irem;
    rate    = r.reactivation(ib)*(1-p.TPTeff);
    mLTBI(destin, source) = mLTBI(destin, source) + rate;

    % --- Initiation of treatment
    sources = [Irec, Irem];
    destin  = Tx;
    rate    = r.gamma;
    m(destin, sources) = m(destin, sources) + rate;
    
    % --- Treatment outcomes 
    source  = Tx;
    destins =  [Rlo,    Rhi];
    rates   = [r.Tx, r.ltfu];
    m(destins, source) = m(destins, source) + rates';

    % --- Self-cure
    sources = intersect(s.infectious, s.(born));
    destin  = Rhi;
    rate    = r.selfcure;
    m(destin, sources) = m(destin, sources) + rate;

    % --- Relapse
    sources = [Rlo, Rhi, R];
    destin  = Irem;
    rates   = r.relapse;
    m(destin, sources) = m(destin, sources) + rates;

    % --- Stabilisation of relapse risk
    sources = [Rlo Rhi];
    destin  = R;
    rates   = 0.5;
    m(destin, sources) = m(destin, sources) + rates;

    % --- Initiation of TPT
    source = Lf;
    destin = Pf;
    rate   = r.TPT(ib);
    m(destin, source) = m(destin, source) + rate;

    source = Lf_imp;
    destin = Pf_imp;
    rate   = r.TPT(ib);
    m(destin, source) = m(destin, source) + rate;

    source = Ls;
    destin = Ps;
    rate   = r.TPT(ib);
    m(destin, source) = m(destin, source) + rate;

    % --- Case-finding
    sources = [Irec, Irem];
    destin  = Tx;
    rate    = r.ACF(ib);
    m(destin, sources) = m(destin, sources) + rate;

end

% Transition from recent to long-term migrant status (over 5-year period)
sources = s.migr_rect;
destins = s.migr_long;
inds = sub2ind([i.nstates, i.nstates], destins, sources);
m(inds) = m(inds) + 1/5;

M.lin = sparse(m - diag(sum(m,1)));
M.LTBIlin = sparse(mLTBI - diag(sum(mLTBI,1)));


% --- Nonlinear component -------------------------------------------------
for ib = 1:length(gps.born)
    born = gps.born{ib};

    m = zeros(i.nstates);
    susinds = intersect([s.U, s.Lf, s.Lf_imp, s.Ls, s.Rlo, s.Rhi, s.R], s.(born));
    m(i.Lf.(born), susinds) = 1;

    imminds = [s.Lf, s.Ls, s.Rlo, s.Rhi, s.R];
    m(:,imminds) = m(:,imminds)*(1-p.imm);

    % Additionally, reinfection amongst people on TPT
    susinds = intersect([s.Pf_imp, s.Ps], s.(born));
    m(i.Pf.(born), susinds) = (1-p.imm);

    M.nlin.(born) = sparse(m - diag(sum(m,1)));     
end


% --- Fractions for different migrant entry states ------------------------

getinds = @(st) intersect(s.migr_rect, s.(st));

prev_in_migr = 0.003;

m = zeros(i.nstates,1);

m(i.U.migr_rect)     = (1-p.LTBIinmigr-prev_in_migr);
m(getinds('Lf_imp')) = p.LTBIinmigr*(1-p.migrTPT)*p.migr_Lf;
m(getinds('Ls'))     = p.LTBIinmigr*(1-p.migrTPT)*(1-p.migr_Lf);
m(getinds('Pf_imp')) = p.LTBIinmigr*p.migrTPT*p.migr_Lf;
m(getinds('Ps'))     = p.LTBIinmigr*p.migrTPT*(1-p.migr_Lf);
m(getinds('Irem'))   = prev_in_migr;

M.migrentries = sparse(m);


% --- Force of infection --------------------------------------------------

getinds = @(st) intersect(s.infectious, s.(st));

m = zeros(2,i.nstates);                                                
m(1,getinds('dom'))  = contmat(1,1);
m(1,getinds('migr')) = contmat(1,2);
m(2,getinds('dom'))  = contmat(2,1);
m(2,getinds('migr')) = contmat(2,2);

% Include infectiousness
m = m*r.beta;
M.lam = sparse(m);


% --- Mortality -----------------------------------------------------------
m = zeros(i.nstates,2);
m(:,1)            = 1/83;
m(s.migr,1)       = m(s.migr,1) + r.migrout;
m(s.infectious,2) = r.muTB;
m(s.Tx,2)         = r.muTx;
M.mort            = sparse(m);
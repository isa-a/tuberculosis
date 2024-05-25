function M = make_model(p, r, i, s, gps)

m = zeros(i.nstates);
for istr = 1:length(gps.strains)
    strain = gps.strains{istr};
    
    % --- Get the linear rates ------------------------------------------------
    L    = i.L.(strain);
    I    = i.I.(strain);
    Dxs  = i.Dx.(strain);
    Txs  = i.Tx.(strain);
    Tx2s = i.Tx2.(strain);
    E    = i.E.(strain);
    Rlo  = i.Rlo.(strain);
    Rhi  = i.Rhi.(strain);
    R    = i.R.(strain);
    
    % --- Reactivation
    source = L; destin = I; rate = r.reactivation;
    m(destin, source) = m(destin, source) + rate;
    
    % --- Primary careseeking, including access to public sector care
    source = I; destin = Dxs.pu; rate = r.careseeking*p.pu + r.access;
    m(destin, source) = m(destin, source) + rate;
    
    source = I; destin = Dxs.pr; rate = r.careseeking*(1-p.pu);
    m(destin, source) = m(destin, source) + rate;
    
    for ip = 1:length(gps.sectors)
        prov = gps.sectors{ip};
        Dx = Dxs.(prov); Tx = Txs.(prov); Tx2 = Tx2s.(prov);
        
        % --- Diagnosis
        ismdr   = strcmp(strain, 'MDR');
        pFLinit = p.Dx(ip)*p.Tx_init(ip)*(1 - ismdr*p.MDR_rec(ip));
        pSLinit = p.Dx(ip)*p.Tx_init(ip)*ismdr*p.MDR_rec(ip);
        p_ltfu  = 1-p.Dx(ip)*p.Tx_init(ip);
        
        source  = Dx;
        destins =          [Tx,       Tx2,      E];
        rates   = r.Dx*[pFLinit,  pSLinit,  p_ltfu];
        m(destins, source) = m(destins, source) + rates';
        
        % --- FL Treatment
        pFLcure = p.cure(ip)*(1-ismdr);
        pSLtran = p.SL_trans(ip)*ismdr;
        rMDRacq = r.MDR_acqu*(1-ismdr);
        
        source  = Tx;
        destins = [Rlo            Tx2                        E                              Rhi,             i.Tx.MDR.(prov)];
        rates   = [r.Tx*pFLcure,  r.Tx*(1-pFLcure)*pSLtran,  r.Tx*(1-pFLcure)*(1-pSLtran),  r.default(ip),   rMDRacq];
        m(destins, source) = m(destins, source) + rates';        
        
        % --- SL Treatment
        source  = Tx2;
        destins = [Rlo                 E                       Rhi];
        rates   = [r.Tx2*p.cure2(ip),  r.Tx2*(1-p.cure2(ip)),  r.default2(ip)];
        m(destins, source) = m(destins, source) + rates';
        
    end
    
    % --- Secondary careseeking
    source = E; destin = Dxs.pu; rate = r.careseeking2*p.pu;
    m(destin, source) = m(destin, source) + rate;
    
    source = E; destin = Dxs.pr; rate = r.careseeking2*(1-p.pu);
    m(destin, source) = m(destin, source) + rate;
    
    % --- Relapse
    sources = [Rlo Rhi R];
    destin  = I;
    rates   = r.relapse;
    m(destin, sources) = m(destin, sources) + rates;
    
    sources = [Rlo Rhi];
    destin  = R;
    rates   = 0.5;
    m(destin, sources) = m(destin, sources) + rates;
    
    % --- Self cure
    sources = intersect(s.infectious,s.(strain));
    destin  = Rhi;
    rates   = r.self_cure;
    m(destin, sources) = m(destin, sources) + rates;

end

M.lin = sparse(m - diag(sum(m,1)));


% --- Get the nonlinear rates ---------------------------------------------

% --- Allocating transitions
for istr = 1:length(gps.strains)
    strain = gps.strains{istr};
    
    m = zeros(i.nstates);
    U = i.U; L = i.L.(strain); I = i.I.(strain);
    
    m(I, [U, s.L, s.Rlo, s.Rhi, s.R]) = p.Fast;
    m(L, [U, s.L, s.Rlo, s.Rhi, s.R]) = 1-p.Fast;
    
    m(:,[s.L, s.Rlo, s.Rhi, s.R]) = m(:,[s.L, s.Rlo, s.Rhi, s.R])*p.imm;
    M.nlin.(strain) = sparse(m - diag(sum(m,1)));
end

% --- Getting force-of-infection
m = zeros(2,i.nstates);
m(1,intersect(s.infectious,s.DS))  = r.beta(1);
m(2,intersect(s.infectious,s.MDR)) = r.beta(2);
m(:,setdiff(s.infectious,s.I)) = m(:,setdiff(s.infectious,s.I))*p.kappa;

M.lambda = sparse(m);


% --- Get the mortality rates
m = r.mort*ones(1,i.nstates);
m([s.infectious, intersect(s.Tx, s.pr)]) = r.mort_TB;
M.mortvec = m';
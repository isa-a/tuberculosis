function M = make_model(p,r,i,s,gps,contmat)

m = zeros(i.nstates);
m2 = zeros(i.nstates);

for ia = 1:length(gps.age)
    age = gps.age{ia};
    
    for ib = 1:length(gps.born)
        born = gps.born{ib};

        for ih = 1:length(gps.hiv)
            hiv = gps.hiv{ih};
            geti = @(st) i.(st).(age).(born).(hiv);
        
            Lf  = geti('Lf');
            Ls  = geti('Ls');
            Pf  = geti('Pf');
            Ps  = geti('Ps');
            I   = geti('I');
            I2  = geti('I2');
            Tx  = geti('Tx');
            Tx2 = geti('Tx2');
            Rlo = geti('Rlo');
            Rhi = geti('Rhi');
            R   = geti('R');
            
            % Progression from 'fast' latent
            source  = Lf;
            destin  = I;
            % rate    = r.progression(ib);
            rate    = r.progression(ia, ib, ih);
            m(destin, source) = m(destin, source) + rate;
    
            source  = Pf;
            destin  = I2;
            %rate    = r.progression(ib)*(1-p.TPTeff(is));
            rate    = r.progression(ia, ib, ih)*(1-p.TPTeff);
            m(destin, source) = m(destin, source) + rate;
    
            % Stabilisation of 'fast' to 'slow' latent
            source = Lf;
            destin = Ls;
            rate   = r.LTBI_stabil;
            m(destin, source) = m(destin, source) + rate;
    
            source = Pf;
            destin = Ps;
            rate   = r.LTBI_stabil;
            m(destin, source) = m(destin, source) + rate;
    
            % Reactivation of 'slow' latent
            source  = Ls;
            destin  = I;
            %rate    = r.reactivation(ib);
            rate    = r.reactivation(ia, ib, ih);
            m(destin, source) = m(destin, source) + rate;
    
            source  = Ps;
            destin  = I;
            %rate    = r.reactivation(ib)*(1-p.TPTeff(is));
            rate    = r.reactivation(ia, ib, ih)*(1-p.TPTeff);
            m(destin, source) = m(destin, source) + rate;
    
            % Initiation of treatment
            source  = I;
            destins = [Tx,                 Tx2,         Rhi];
            rates   = [r.gamma(ia), r.gamma(ia), r.self_cure];
            m(destins, source) = m(destins, source) + rates';
    
            source  = I2;
            destins = [Tx,                 Tx2,         Rhi];
            rates   = [r.gamma(ia), r.gamma(ia), r.self_cure];
            m(destins, source) = m(destins, source) + rates';
    
            % Treatment completion or interruption
            source  = Tx;
            destins = [Rlo Rhi];
            rates   = [r.Tx, r.ltfu];
            m(destins, source) = m(destins, source) + rates';

    
            % Second-line treatment
            source  = Tx2;
            destins = [Rlo Rhi];
            rates   = [r.Tx2, r.ltfu2];
            m(destins, source) = m(destins, source) + rates';
    
            % Relapse
            sources = [Rlo Rhi R];
            destin  = I2;
            rates   = r.relapse;
            m(destin, sources) = m(destin, sources) + rates;
    
            % Stabilisation of relapse risk
            sources = [Rlo Rhi];
            destin  = R;
            rates   = 0.5;
            m(destin, sources) = m(destin, sources) + rates;
    
            % Initiation of TPT
            source = Lf;
            destin = Pf;
            rate   = r.TPT(ib);
            m(destin, source) = m(destin, source) + rate;
    
            source = Ls;
            destin = Ps;
            rate   = r.TPT(ib);
            m(destin, source) = m(destin, source) + rate;
    
            % Case-finding
            sources = [I I2];
            destin  = Tx;
            rate    = r.ACF(ib);
            m(destin, sources) = m(destin, sources) + rate;
    
            sources = [I I2];
            destin  = Tx2;
            rate    = r.ACF(ib);
            m(destin, sources) = m(destin, sources) + rate;
    
    
            source = I2;
            destin = Tx;
            rate   = r.ACF2(ib);
            m(destin, source) = m(destin, source) + rate;
    
            source = I2;
            destin = Tx2;
            rate   = r.ACF2(ib);
            m(destin, source) = m(destin, source) + rate;
        end
    end
end

% Transition from recent to long-term migrant status (over 5-year period)
% sources = s.migr_rect;
% destins = s.migr_long;
% inds = sub2ind([i.nstates, i.nstates], destins, sources);
% m(inds) = m(inds) + 1/5;

% Transition from dom to vulnerable population
% sources = s.dom;
% destins = s.vuln;
% inds = sub2ind([i.nstates, i.nstates], destins, sources);
% m(inds) = m(inds) + r.vuln;

% --- Ageing process
% sources = s.ch;
% destins = s.ad;
% inds = sub2ind([i.nstates, i.nstates], destins, sources);
% m(inds) = m(inds) + r.ageing;

% --- HIV acquisition
sources = s.neg;
destins = s.pos;
inds    = sub2ind(size(m), destins, sources);
rates   = 1;
m2(inds) = m2(inds) + rates;

% --- ART initiation
sources = s.pos;
destins = s.art;
inds    = sub2ind(size(m), destins, sources);
rates   = r.ART_init;
m(inds) = m(inds) + rates;

% --- Combine
M.lin    = sparse(m - diag(sum(m,1)));
M.linHIV = sparse(m2 - diag(sum(m2,1)));


% --- Nonlinear component -------------------------------------------------

for ia = 1:length(gps.age)
    age = gps.age{ia};
    m = zeros(i.nstates); % new added here instead of after hiv
    for ib = 1:length(gps.born)
        born = gps.born{ib};
        for ih = 1:length(gps.hiv)
            hiv = gps.hiv{ih};

            susinds = intersect(intersect(intersect([s.U, s.Lf, s.Ls, s.Rlo, s.Rhi, s.R],s.(age)),s.(born)), s.(hiv));
            m(i.Lf.(age).(born).(hiv), susinds) = 1;
            
            imminds = [s.Lf, s.Ls, s.Rlo, s.Rhi, s.R];
            m(:,imminds) = m(:,imminds)*(1-p.imm);
            
            M.nlin.(age) = sparse(m - diag(sum(m,1)));     % <--- Make sure all of these are used in goveqs_basis, multiplied by relevant elements of lambda
        end
    end
end


% --- Force of infection --------------------------------------------------

% getinds = @(st1, st2) intersect(intersect(s.infectious, s.(st1)), s.(st2));
% contmat(end,end) = contmat(end,end);
% 
% m = zeros(2,i.nstates);                                                     % Rows: 1.Dom DS 2.Dom RR 3.Migr DS 4.Migr RR 5.Vuln DS 6.Vuln RR
%                                                                             % no RR Rows: 1.Dom DS 2.Migr DS 3.Vuln DS   

% m(1,getinds('ch', 'dom')) = contmat(1,1);
% m(1,getinds('ad', 'dom')) = contmat(1,2);
% m(1,getinds('ad', 'dom', 'ds')) = contmat(1,2);                                  
% m(1,getinds('ch', 'vuln','ds')) = contmat(1,3);
% m(1,getinds('ad', 'vuln','ds')) = contmat(1,4);

% m(2,getinds('ch', 'dom')) = contmat(2,1);
% m(2,getinds('ad', 'dom')) = contmat(2,2);
% m(2,getinds('ad', 'dom','ds')) = contmat(2,2);
% m(2,getinds('ch', 'vuln','ds')) = contmat(2,3);
% m(2,getinds('ad', 'vuln','ds')) = contmat(2,4);

% m(3,getinds('ch', 'dom', 'ds')) = contmat(3,1);
% m(3,getinds('ad', 'dom','ds')) = contmat(3,2);
% m(3,getinds('ch', 'vuln','ds')) = contmat(3,3);
% m(3,getinds('ad', 'vuln','ds')) = contmat(3,4);
% 
% m(4,getinds('ch', 'dom', 'ds')) = contmat(4,1);
% m(4,getinds('ad', 'dom','ds')) = contmat(4,2);
% m(4,getinds('ch', 'vuln','ds')) = contmat(4,3);
% m(4,getinds('ad', 'vuln','ds')) = contmat(4,4);


% Include infectiousness
m = m*r.beta;
M.lam = sparse(m);

% Additional matrix to help keep track of numbers in each group 
m = zeros(2,i.nstates);
% m(1, intersect(s.ch, s.dom))  = 1;  
m(1, intersect(s.ad, s.dom))  = 1; 
% m(3, intersect(s.ch, s.migr)) = 1; 
% m(4, intersect(s.ad, s.migr)) = 1; 
% m(3, intersect(s.ch, s.vuln)) = 1; 
% m(4, intersect(s.ad, s.vuln)) = 1;  
M.denvec = sparse(m);


% --- Mortality -----------------------------------------------------------
m = zeros(i.nstates,2);
% m(s.ch,1)         = 0;
m(s.ad,1)         = 1/72;
m(s.pos,1)        = m(s.pos,1) + r.HIV_mort;
%m(:,1)            = 1/83;
% m(s.vuln,1)       = 1/55;
m(s.infectious,2) = r.muTB;

M.mort            = sparse(m);

% --- Mortality -----------------------------------------------------------

% m = zeros(i.nstates,1);
% m(s.for) = r.migr;
% m(intersect(s.Lf, s.for)) = r.migr*p.kLf;
% M.migration = sparse(m);
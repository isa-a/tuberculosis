function M = make_model(p,r,i,s,gps,contmat)

m = zeros(i.nstates);


for ia = 1:length(gps.age)
    age = gps.age{ia};

    for is = 1:length(gps.strains)
        strain = gps.strains{is};
        ismdr = strcmp(strain,'rr');
    
        for ib = 1:length(gps.born)
            born = gps.born{ib};

            for ih = 1:length(gps.hiv)
                hiv = gps.hiv{ih};
                geti = @(st) i.(st).(age).(born).(strain).(hiv);
                
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
                rate    = r.progression(ia, ib, ih)*(1-p.TPTeff(is));
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
                rate    = r.reactivation(ia, ib, ih)*(1-p.TPTeff(is));
                m(destin, source) = m(destin, source) + rate;
        
                % Initiation of treatment
                pSLinit = ismdr*p.RRrec;
                source  = I;
                destins =                      [Tx,                 Tx2,         Rhi];
                rates   = [r.gamma(ia)*(1-pSLinit), r.gamma(ia)*pSLinit, r.self_cure];
                m(destins, source) = m(destins, source) + rates';
        
                source  = I2;
                destins =                      [Tx,                 Tx2,         Rhi];
                rates   = [r.gamma(ia)*(1-pSLinit), r.gamma(ia)*pSLinit, r.self_cure];
                m(destins, source) = m(destins, source) + rates';
        
                % Treatment completion or interruption
                source  = Tx;
                destins = [Rlo Rhi];
                rates   = [r.Tx, r.ltfu];
                m(destins, source) = m(destins, source) + rates';
        
                % Acquisition of drug resistance while on first-line treatment
                if ~ismdr
                    source = Tx;
                    destin = i.Tx.(age).(born).rr.(hiv);                                   % <--- Include age stratification
                    rate   = r.RR_acqu;
                    m(destin, source) = m(destin, source) + rate;
                end
        
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
                rate    = r.ACF(ib)*(1-pSLinit);
                m(destin, sources) = m(destin, sources) + rate;
        
                sources = [I I2];
                destin  = Tx2;
                rate    = r.ACF(ib)*pSLinit;
                m(destin, sources) = m(destin, sources) + rate;
        
        
                source = I2;
                destin = Tx;
                rate   = r.ACF2(ib)*(1-pSLinit);
                m(destin, source) = m(destin, source) + rate;
        
                source = I2;
                destin = Tx2;
                rate   = r.ACF2(ib)*pSLinit;
                m(destin, source) = m(destin, source) + rate;
            end
        end
    end
end

% Transition from recent to long-term migrant status (over 5-year period)
sources = s.migr_rect;
destins = s.migr_long;
inds = sub2ind([i.nstates, i.nstates], destins, sources);
m(inds) = m(inds) + 1/5;

% Transition from dom to vulnerable population
sources = s.dom;
destins = s.vuln;
inds = sub2ind([i.nstates, i.nstates], destins, sources);
m(inds) = m(inds) + r.vuln;

% --- Ageing process
sources = s.ch;
destins = s.ad;
inds = sub2ind([i.nstates, i.nstates], destins, sources);
m(inds) = m(inds) + r.ageing;

M.lin = sparse(m - diag(sum(m,1)));


% --- Nonlinear component -------------------------------------------------

for ia = 1:length(gps.age)
    age = gps.age{ia};
    for is = 1:length(gps.strains)
        strain = gps.strains{is};
        for ib = 1:length(gps.born)
            born = gps.born{ib};
    
            m = zeros(i.nstates);
            susinds = intersect(intersect([s.U, s.Lf, s.Ls, s.Rlo, s.Rhi, s.R],s.(age)),s.(born));
            m(i.Lf.(age).(born).(strain), susinds) = 1;
            
            imminds = [s.Lf, s.Ls, s.Rlo, s.Rhi, s.R];
            m(:,imminds) = m(:,imminds)*(1-p.imm);
            
            M.nlin.(age).(born).(strain) = sparse(m - diag(sum(m,1)));     % <--- Make sure all of these are used in goveqs_basis, multiplied by relevant elements of lambda
        end
    end
end

% --- Fractions for different migrant entry states ------------------------

getindsch = @(st1, st2) intersect(intersect(intersect(s.migr_rect, s.(st1)), s.(st2)),s.ch);
getindsad = @(st1, st2) intersect(intersect(intersect(s.migr_rect, s.(st1)), s.(st2)),s.ad);
m = zeros(i.nstates,1);

m(i.U.ch.migr_rect) = (1-p.LTBI_in_migrch)*p.ch_in_migr;
m(i.U.ad.migr_rect) = (1-p.LTBI_in_migrad)*(1-p.ch_in_migr);

m(getindsch('Lf','ds')) = p.LTBI_in_migrch*p.ch_in_migr*(1-p.migrTPT)*0.02;

m(getindsch('Ls','ds')) = p.LTBI_in_migrch*p.ch_in_migr*(1-p.migrTPT)*0.98;

m(getindsch('Pf','ds')) = p.LTBI_in_migrch*p.ch_in_migr*p.migrTPT*0.02;

m(getindsch('Ps','ds')) = p.LTBI_in_migrch*p.ch_in_migr*p.migrTPT*0.98;


m(getindsad('Lf','ds')) = p.LTBI_in_migrad*(1-p.ch_in_migr)*(1-p.migrTPT)*0.02;

m(getindsad('Ls','ds')) = p.LTBI_in_migrad*(1-p.ch_in_migr)*(1-p.migrTPT)*0.98;

m(getindsad('Pf','ds')) = p.LTBI_in_migrad*(1-p.ch_in_migr)*p.migrTPT*0.02;

m(getindsad('Ps','ds')) = p.LTBI_in_migrad*(1-p.ch_in_migr)*p.migrTPT*0.98;


% m(getindsch('Lf','ds')) = p.LTBI_in_migrch*p.ch_in_migr*(1-p.migrTPT)*(1-p.RR_in_migr)*0.02;
% m(getindsch('Lf','rr')) = p.LTBI_in_migrch*p.ch_in_migr*(1-p.migrTPT)*p.RR_in_migr*0.02;
% m(getindsch('Ls','ds')) = p.LTBI_in_migrch*p.ch_in_migr*(1-p.migrTPT)*(1-p.RR_in_migr)*0.98;
% m(getindsch('Ls','rr')) = p.LTBI_in_migrch*p.ch_in_migr*(1-p.migrTPT)*p.RR_in_migr*0.98;
% m(getindsch('Pf','ds')) = p.LTBI_in_migrch*p.ch_in_migr*p.migrTPT*(1-p.RR_in_migr)*0.02;
% m(getindsch('Pf','rr')) = p.LTBI_in_migrch*p.ch_in_migr*p.migrTPT*p.RR_in_migr*0.02;
% m(getindsch('Ps','ds')) = p.LTBI_in_migrch*p.ch_in_migr*p.migrTPT*(1-p.RR_in_migr)*0.98;
% m(getindsch('Ps','rr')) = p.LTBI_in_migrch*p.ch_in_migr*p.migrTPT*p.RR_in_migr*0.98;
% 
% m(getindsad('Lf','ds')) = p.LTBI_in_migrad*(1-p.ch_in_migr)*(1-p.migrTPT)*(1-p.RR_in_migr)*0.02;
% m(getindsad('Lf','rr')) = p.LTBI_in_migrad*(1-p.ch_in_migr)*(1-p.migrTPT)*p.RR_in_migr*0.02;
% m(getindsad('Ls','ds')) = p.LTBI_in_migrad*(1-p.ch_in_migr)*(1-p.migrTPT)*(1-p.RR_in_migr)*0.98;
% m(getindsad('Ls','rr')) = p.LTBI_in_migrad*(1-p.ch_in_migr)*(1-p.migrTPT)*p.RR_in_migr*0.98;
% m(getindsad('Pf','ds')) = p.LTBI_in_migrad*(1-p.ch_in_migr)*p.migrTPT*(1-p.RR_in_migr)*0.02;
% m(getindsad('Pf','rr')) = p.LTBI_in_migrad*(1-p.ch_in_migr)*p.migrTPT*p.RR_in_migr*0.02;
% m(getindsad('Ps','ds')) = p.LTBI_in_migrad*(1-p.ch_in_migr)*p.migrTPT*(1-p.RR_in_migr)*0.98;
% m(getindsad('Ps','rr')) = p.LTBI_in_migrad*(1-p.ch_in_migr)*p.migrTPT*p.RR_in_migr*0.98;

M.migrentries = sparse(m);


% --- Force of infection --------------------------------------------------

getinds = @(st1, st2, st3) intersect(intersect(intersect(s.infectious, s.(st1)), s.(st2)), s.(st3));
contmat(end,end) = contmat(end,end);

m = zeros(6,i.nstates);                                                   % Rows: 1.Dom DS 2.Dom RR 3.Migr DS 4.Migr RR 5.Vuln DS 6.Vuln RR
                                                                            % no RR Rows: 1.Dom DS 2.Migr DS 3.Vuln DS   
% for ii = 1:6
%     m(ii,getinds('ch', 'dom', 'ds')) = contmat(ii,1);                              % no vuln Rows: 1.Dom DS 2.Migr DS
%     m(ii,getinds('ad', 'dom', 'ds')) = contmat(ii,2);
%     m(ii,getinds('ch', 'migr','ds')) = contmat(ii,3);
%     m(ii,getinds('ad', 'migr','ds')) = contmat(ii,4);
%     m(ii,getinds('ch', 'vuln','ds')) = contmat(ii,5);
%     m(ii,getinds('ad', 'vuln','ds')) = contmat(ii,6);
% end

m(1,getinds('ch', 'dom', 'ds')) = contmat(1,1);                              % no vuln Rows: 1.Dom DS 2.Migr DS 
m(1,getinds('ad', 'dom', 'ds')) = contmat(1,2);                                  
m(1,getinds('ch', 'migr','ds')) = contmat(1,3);
m(1,getinds('ad', 'migr','ds')) = contmat(1,4);
m(1,getinds('ch', 'vuln','ds')) = contmat(1,5);
m(1,getinds('ad', 'vuln','ds')) = contmat(1,6);

% m(2,getinds('ch', 'dom', 'rr')) = contmat(1,1);
% m(2,getinds('ad', 'dom','rr')) = contmat(1,2);
% m(2,getinds('ch', 'migr','rr')) = contmat(1,3);
% m(2,getinds('ad', 'migr','rr')) = contmat(1,4);
% m(2,getinds('ch', 'vuln','rr')) = contmat(1,5);
% m(2,getinds('ad', 'vuln','rr')) = contmat(1,6);

m(2,getinds('ch', 'dom', 'ds')) = contmat(2,1);
m(2,getinds('ad', 'dom','ds')) = contmat(2,2);
m(2,getinds('ch', 'migr','ds')) = contmat(2,3);
m(2,getinds('ad', 'migr','ds')) = contmat(2,4);
m(2,getinds('ch', 'vuln','ds')) = contmat(2,5);
m(2,getinds('ad', 'vuln','ds')) = contmat(2,6);

% m(4,getinds('ch', 'dom', 'rr')) = contmat(2,1);
% m(4,getinds('ad', 'dom','rr')) = contmat(2,2);
% m(4,getinds('ch', 'migr','rr')) = contmat(2,3);
% m(4,getinds('ad', 'migr','rr')) = contmat(2,4);
% m(4,getinds('ch', 'vuln','rr')) = contmat(2,5);
% m(4,getinds('ad', 'vuln','rr')) = contmat(2,6);

m(3,getinds('ch', 'dom', 'ds')) = contmat(3,1);
m(3,getinds('ad', 'dom','ds')) = contmat(3,2);
m(3,getinds('ch', 'migr','ds')) = contmat(3,3);
m(3,getinds('ad', 'migr','ds')) = contmat(3,4);
m(3,getinds('ch', 'vuln','ds')) = contmat(3,5);
m(3,getinds('ad', 'vuln','ds')) = contmat(3,6);

% m(6,getinds('ch', 'dom', 'rr')) = contmat(3,1);
% m(6,getinds('ad', 'dom','rr')) = contmat(3,2);
% m(6,getinds('ch', 'migr','rr')) = contmat(3,3);
% m(6,getinds('ad', 'migr','rr')) = contmat(3,4);
% m(6,getinds('ch', 'vuln','rr')) = contmat(3,5);
% m(6,getinds('ad', 'vuln','rr')) = contmat(3,6);

m(4,getinds('ch', 'dom', 'ds')) = contmat(4,1);
m(4,getinds('ad', 'dom','ds')) = contmat(4,2);
m(4,getinds('ch', 'migr','ds')) = contmat(4,3);
m(4,getinds('ad', 'migr','ds')) = contmat(4,4);
m(4,getinds('ch', 'vuln','ds')) = contmat(4,5);
m(4,getinds('ad', 'vuln','ds')) = contmat(4,6);

% m(8,getinds('ch', 'dom', 'rr')) = contmat(4,1);
% m(8,getinds('ad', 'dom','rr')) = contmat(4,2);
% m(8,getinds('ch', 'migr','rr')) = contmat(4,3);
% m(8,getinds('ad', 'migr','rr')) = contmat(4,4);
% m(8,getinds('ch', 'vuln','rr')) = contmat(4,5);
% m(8,getinds('ad', 'vuln','rr')) = contmat(4,6);

m(5,getinds('ch', 'dom', 'ds')) = contmat(5,1);
m(5,getinds('ad', 'dom','ds')) = contmat(5,2);
m(5,getinds('ch', 'migr','ds')) = contmat(5,3);
m(5,getinds('ad', 'migr','ds')) = contmat(5,4);
m(5,getinds('ch', 'vuln','ds')) = contmat(5,5);
m(5,getinds('ad', 'vuln','ds')) = contmat(5,6);

% m(10,getinds('ch', 'dom', 'rr')) = contmat(5,1);
% m(10,getinds('ad', 'dom','rr')) = contmat(5,2);
% m(10,getinds('ch', 'migr','rr')) = contmat(5,3);
% m(10,getinds('ad', 'migr','rr')) = contmat(5,4);
% m(10,getinds('ch', 'vuln','rr')) = contmat(5,5);
% m(10,getinds('ad', 'vuln','rr')) = contmat(5,6);

m(6,getinds('ch', 'dom', 'ds')) = contmat(6,1);
m(6,getinds('ad', 'dom','ds')) = contmat(6,2);
m(6,getinds('ch', 'migr','ds')) = contmat(6,3);
m(6,getinds('ad', 'migr','ds')) = contmat(6,4);
m(6,getinds('ch', 'vuln','ds')) = contmat(6,5);
m(6,getinds('ad', 'vuln','ds')) = contmat(6,6);

% m(12,getinds('ch', 'dom', 'rr')) = contmat(6,1);
% m(12,getinds('ad', 'dom','rr')) = contmat(6,2);
% m(12,getinds('ch', 'migr','rr')) = contmat(6,3);
% m(12,getinds('ad', 'migr','rr')) = contmat(6,4);
% m(12,getinds('ch', 'vuln','rr')) = contmat(6,5);
% m(12,getinds('ad', 'vuln','rr')) = contmat(6,6);

% Include infectiousness
m = m*r.beta;
% Discount overall infectiousness of RR
%m(:,intersect(s.rr,s.infectious)) = m(:,intersect(s.rr,s.infectious))*p.relbeta_RR;

M.lam = sparse(m);

% Additional matrix to help keep track of numbers in each group 
m = zeros(6,i.nstates);
m(1, intersect(s.ch, s.dom))  = 1;  
m(2, intersect(s.ad, s.dom))  = 1; 
m(3, intersect(s.ch, s.migr)) = 1; 
m(4, intersect(s.ad, s.migr)) = 1; 
m(5, intersect(s.ch, s.vuln)) = 1; 
m(6, intersect(s.ad, s.vuln)) = 1;  
M.denvec = sparse(m);


% --- Mortality -----------------------------------------------------------
m = zeros(i.nstates,2);
m(s.ch,1)         = 0;
m(s.ad,1)         = 1/72;
%m(:,1)            = 1/83;
m(s.vuln,1)       = 1/55;
m(s.infectious,2) = r.muTB;
M.mort            = sparse(m);

% --- Mortality -----------------------------------------------------------

% m = zeros(i.nstates,1);
% m(s.for) = r.migr;
% m(intersect(s.Lf, s.for)) = r.migr*p.kLf;
% M.migration = sparse(m);
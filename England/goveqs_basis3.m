function out = goveqs_basis2(t, in, i, s, M, rin_vec, agg, sel, r, p, equil)

out = zeros(length(in),1);
invec = in(1:i.nstates);

% Prepare population denominators
tmp = M.denvec*invec;
den = sum(tmp.*M.denvec,1)';
den(den==0) = Inf;

% New infections
% lam = M.lam*(invec./den)*(1-p.betadec)^(max((t-2010),0));
lam = M.lam*(invec./den);
% lam = M.lam*(invec./sum(invec))*(1-p.betadec)^(max((t-2010),0));

% Indices for lambda entries
% 1. ch dom  ds
% 2. ch dom  rr
% 3. ad dom  ds
% 4. ad dom  rr
% 5. ch migr ds
% 6. ch migr rr
% 7. ad migr ds
% 8. ad migr rr
% 9. ch vuln ds
% 10.ch vuln rr
% 11.ad vuln ds
% 12.ad vuln rr

% Indices for NO RR lambda entries
% 1. ch dom  ds
% 2. ad dom  ds
% 3. ch migr ds
% 4. ad migr ds
% 5. ch vuln ds
% 6. ad vuln ds


% Indices for NO RR NO VULN lambda entries
% 1. ch dom  ds
% 2. ad dom  ds
% 3. ch migr ds
% 4. ad migr ds


% allmat = M.lin + ...
%         lam(1)*M.nlin.ch.dom.ds           + lam(2)*M.nlin.ad.dom.ds + ...
%         lam(3)*M.nlin.ch.migr_rect.ds     + lam(4)*M.nlin.ad.migr_rect.ds + ...
%         lam(3)*M.nlin.ch.migr_long.ds     + lam(4)*M.nlin.ad.migr_long.ds;

allmat = M.lin + M.LTBIlin*(1-p.LTBIdec)^min(max((t-2010),0),7) + ...
        lam(1)*M.nlin.ad.dom.ds + ...
        lam(2)*M.nlin.ad.migr_rect.ds + ...
        lam(2)*M.nlin.ad.migr_long.ds; % + ...
        % lam(5)*M.nlin.ch.vuln.ds          + lam(6)*M.nlin.ad.vuln.ds;


% Full model
% allmat = M.lin + ...
%         lam(1)*M.nlin.ch.dom.ds           + lam(2)*M.nlin.ch.dom.rr + ...
%         lam(3)*M.nlin.ad.dom.ds           + lam(4)*M.nlin.ad.dom.rr + ...
%         lam(5)*M.nlin.ch.migr_rect.ds     + lam(6)*M.nlin.ch.migr_rect.rr + ...
%         lam(5)*M.nlin.ch.migr_long.ds     + lam(6)*M.nlin.ch.migr_long.rr + ...
%         lam(7)*M.nlin.ad.migr_rect.ds     + lam(8)*M.nlin.ad.migr_rect.rr + ...
%         lam(7)*M.nlin.ad.migr_long.ds     + lam(8)*M.nlin.ad.migr_long.rr + ...
%         lam(9)*M.nlin.ch.vuln.ds          + lam(10)*M.nlin.ch.vuln.rr + ...
%         lam(11)*M.nlin.ad.vuln.ds         + lam(12)*M.nlin.ad.vuln.rr;

out(1:i.nstates) = allmat*invec;

% Mortality
morts = M.mort.*invec;
out(1:i.nstates) = out(1:i.nstates) - sum(morts,2);

% Births into UK population
dom_morts = sum(sum(morts([s.dom],:)));
out(i.U.ad.dom) = out(i.U.ad.dom) + dom_morts;

% method 1 ----------------------------------------------------------------

% years = [2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023];
% vals = [797, 772, 752, 825, 886, 682, 773, 1294, 1316];
%mgrtn_ft = interp1(years, vals ./ vals(years==2015), t, 'linear', 'extrap');

% if equil
%     r_in = rin_vec(1);
% else
%     r_in = rin_vec(min(max(ceil(t-2019),1),length(rin_vec)));
% end
% 
% out(1:i.nstates) = out(1:i.nstates) + r_in*M.migrentries; 

if equil
    r_in = rin_vec(1);
else
    r_in = rin_vec(min(max(round(t) - 2019 + 1,1),length(rin_vec)));
end

out(1:i.nstates) = out(1:i.nstates) + r_in*M.migrentries;



% if t < 2015
%     mgrtn_ft = 1;
%     % Migration out of UK
%     outmigr = r.migr * invec(s.migr) / sum(invec(s.migr));
%     out(s.migr) = out(s.migr) - outmigr;
% 
%     % Migration into UK 
%     migrmorts = sum(sum(morts(s.migr,:)));
%     totalout = sum(outmigr) + migrmorts;
%     out(1:i.nstates) = out(1:i.nstates) + mgrtn_ft * totalout .* M.migrentries;
% else
%     mgrtn_ft = interp1(years, vals ./ vals(years==2015), t, 'linear', 'extrap');
%     % Migration out of UK
%     outmigr = 1 / r.migr;
%     out(s.migr) = out(s.migr) - outmigr;
% 
%     % Migration into UK
%     out(1:i.nstates) = out(1:i.nstates) + mgrtn_ft * sum(outmigr) .* M.migrentries;
% end

% % Migration out of UK
% outmigr = r.migr * invec(s.migr) / sum(invec(s.migr));
% out(s.migr) = out(s.migr) - outmigr;
% 
% % Migration into UK 
% migrmorts = sum(sum(morts(s.migr,:)));
% totalout = sum(outmigr) + migrmorts;
% out(1:i.nstates) = out(1:i.nstates) + mgrtn_ft * totalout .* M.migrentries;

% keyboard;

% %disp('in');
% %disp(sum(out));

% method 2 ----------------------------------------------------------------

% % Migration out of UK
% outmigr = r.migr*invec(s.migr)/sum(invec(s.migr));
% out(s.migr) = out(s.migr) - outmigr;
% 
% % Migration into UK
% migrmorts = sum(sum(morts(s.migr,:)));
% totalout = sum(outmigr) + migrmorts;
% % migration into migr recent states only
% out(s.migr_rect) = out(s.migr_rect) + (totalout / length(s.migr_rect));


if sum(invec)>5
    keyboard;
end


% % Migration
% out(1:i.nstates) = out(1:i.nstates) + M.migration;

% Auxiliaries
out(i.aux.inc)        = agg.inc*(sel.inc.*allmat)*invec;

tmp1                  = agg.incsources*((sel.L2I.*allmat)*invec);
tmp2                  = agg.incsources*((sel.P2I.*allmat)*invec);
tmp3                  = agg.incsources*((sel.R2I.*allmat)*invec);
tmp4                  = agg.incsources*((sel.T2I.*allmat)*invec);
out(i.aux.incsources) = [tmp1; tmp2; tmp3; tmp4];
out(i.aux.mort)       = sum(morts(:,2));
out(i.aux.nTPT)       = sum((sel.nTPT.*allmat)*invec);
out(i.aux.ch_notifs)  = sum((sel.ch_notifs.*allmat)*invec);
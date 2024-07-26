function out = goveqs_basis2(t, in, i, s, M, agg, sel, r, p)

out = zeros(length(in),1);
invec = in(1:i.nstates);

% Prepare population denominators
tmp = M.denvec*invec;
den = sum(tmp.*M.denvec,1)';
den(den==0) = Inf;

% New infections
lam = M.lam*(invec./den)*(1-p.betadec)^(max((t-2010),0));
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

% Full model
allmat = M.lin + ...
        lam(1)*M.nlin.ch.dom.ds           + lam(2)*M.nlin.ch.dom.rr + ...
        lam(3)*M.nlin.ad.dom.ds           + lam(4)*M.nlin.ad.dom.rr + ...
        lam(5)*M.nlin.ch.migr_rect.ds     + lam(6)*M.nlin.ch.migr_rect.rr + ...
        lam(5)*M.nlin.ch.migr_long.ds     + lam(6)*M.nlin.ch.migr_long.rr + ...
        lam(7)*M.nlin.ad.migr_rect.ds     + lam(8)*M.nlin.ad.migr_rect.rr + ...
        lam(7)*M.nlin.ad.migr_long.ds     + lam(8)*M.nlin.ad.migr_long.rr + ...
        lam(9)*M.nlin.ch.vuln.ds          + lam(10)*M.nlin.ch.vuln.rr + ...
        lam(11)*M.nlin.ad.vuln.ds         + lam(12)*M.nlin.ad.vuln.rr;

out(1:i.nstates) = allmat*invec;
%disp('after allmat:');
%disp(sum(out));

% Mortality
morts = M.mort.*invec;
out(1:i.nstates) = out(1:i.nstates) - sum(morts,2);
%disp('after mortality:');
%disp(sum(out));

% Births into UK population
dom_morts = sum(sum(morts([s.dom,s.vuln],:)));
out(i.U.ch.dom) = out(i.U.ch.dom) + dom_morts;
%disp('after births:');
%disp(sum(out));

% Migration out of UK
outmigr = r.migr*invec(s.migr)/sum(invec(s.migr));
outmigr = min(outmigr, invec(s.migr)); % so migration out doesnt exceed migrpop
out(s.migr) = out(s.migr) - outmigr;

% Migration into UK
inmigr = sum(sum(morts(s.migr,:))) + r.migr;
% vec = [1-p.LTBI_in_migr, (1-p.migrTPT)*p.LTBI_in_migr*[0.02 0.98], p.migrTPT*p.LTBI_in_migr*[0.02 0.98]]';
% out(s.migrstates) = out(s.migrstates) + inmigr*vec;
inmigr = min(inmigr, sum(invec(s.migr))); % make sure incoming migration isnt more than number in migr states
out(1:i.nstates) = out(1:i.nstates) + inmigr * M.migrentries;

%disp('after migration in:');
%disp(sum(out));

% final population
% pop = sum(out(1:i.nstates));
% if pop > 5
%     % bring population down by scale of 5
%     scale = 5 / pop;
%     out = out * scale;
% end

% disp('final population sum:');
% disp(sum(out));
end

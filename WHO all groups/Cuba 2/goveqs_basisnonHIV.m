function out = goveqs_basisnonHIV(t, in, i, s, M, agg, sel, r, p)

out = zeros(length(in),1);
invec = in(1:i.nstates);

% Prepare population denominators
tmp = M.denvec*invec;
den = sum(tmp.*M.denvec,1)';
den(den==0) = Inf;

% if t<1980                                                                  
%     rHIV = 0;
% else
%     rHIV = interp1(0:length(r.rHIV)-1, r.rHIV, (t-1980)); 
% end

try
% Normalise by populations
lam = M.lam*(invec./den)*(1-p.betadec)^(max((t-2010),0));
allmat = M.lin + lam(1)*M.nlin.ch + lam(2)*M.nlin.ad;                     
out = allmat*invec;         
catch
   keyboard; 
end


% Mortality and births
morts = M.mort.*invec;
out(1:i.nstates) = out(1:i.nstates) - sum(morts,2);

dom_morts = sum(sum(morts([s.dom,s.vuln],:)));
out(i.U.ch.dom.neg) = out(i.U.ch.dom.neg) + dom_morts;


if sum(invec)>5
    keyboard;
end


% Auxiliaries
out(i.aux.inc)        = agg.inc*(sel.inc.*allmat)*invec;
% tmp1                  = agg.incsources*((sel.L2I.*allmat)*invec);
% tmp2                  = agg.incsources*((sel.P2I.*allmat)*invec);
% tmp3                  = agg.incsources*((sel.R2I.*allmat)*invec);
% tmp4                  = agg.incsources*((sel.T2I.*allmat)*invec);
% out(i.aux.incsources) = [tmp1; tmp2; tmp3; tmp4];
out(i.aux.mort)       = sum(morts(:,2));
out(i.aux.nTPT)       = sum((sel.nTPT.*allmat)*invec);
out(i.aux.ch_notifs)  = sum((sel.ch_notifs.*allmat)*invec);

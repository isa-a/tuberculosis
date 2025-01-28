gps0 = gps;
gps0.vacc  = {'v0'};

[i0, s0, d0, lim] = get_addresses({states1, gps0.vacc, gps0.hiv}, [], [], [], 0);
d0 = char(d0);

s0.infectious  = [s0.Isc, s0.I, s0.E, s0.Dx];
s0.symptomatic = [s0.I, s0.E, s0.Dx];
s0.prevalent   = [s0.infectious, s0.Tx];


% --- Include the auxiliaries ---------------------------------------------
names = {'inc','noti','mort','Dx'};
lgths =     [3,     3,     2,   1];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    i0.aux.(names{ii}) = inds;
    lim = inds(end);
end
i0.nx = lim;


% --- Make aggregators and selectors --------------------------------------

% Selectors for the incidence
tmp = zeros(3,i0.nstates);
tmp(1,intersect(s0.Isc,s0.h0)) = 1;
tmp(2,intersect(s0.Isc,s0.h1)) = 1;
tmp(3,intersect(s0.Isc,s0.hart)) = 1;
agg0.inc = sparse(tmp);

tmp = zeros(i0.nstates);
tmp(s0.Isc,:) = 1;
tmp(s0.h1,s0.h0) = 0;   tmp(s0.h0,s0.h1) = 0; 
tmp(s0.hart,s0.h1) = 0; tmp(s0.h1,s0.hart) = 0;
sel0.inc = tmp - diag(diag(tmp));

% Selectors for notifications
tmp = zeros(3,i0.nstates);
tmp(1,intersect(s0.Tx,s0.h0)) = 1;
tmp(2,intersect(s0.Tx,s0.h1)) = 1;
tmp(3,intersect(s0.Tx,s0.hart)) = 1;
agg0.noti = sparse(tmp);

tmp = zeros(i0.nstates);
tmp(s0.Tx,:) = 1;
tmp(s0.h1,s0.h0) = 0;   tmp(s0.h0,s0.h1) = 0; 
tmp(s0.hart,s0.h1) = 0; tmp(s0.h1,s0.hart) = 0;
sel0.noti = tmp - diag(diag(tmp));

% Selectors for counting presentations to care
tmp = zeros(1,i0.nstates);
tmp(1,s0.Dx) = 1;
agg0.Dx = sparse(tmp);

tmp = zeros(i0.nstates);
tmp(s0.Dx,:) = 1;
tmp(s0.h1,s0.h0) = 0;   tmp(s0.h0,s0.h1) = 0; 
tmp(s0.hart,s0.h1) = 0; tmp(s0.h1,s0.hart) = 0;
sel0.Dx = tmp - diag(diag(tmp));
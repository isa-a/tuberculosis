gps.hiv = {'h0','h1','hart'};
gps.vacc  = {'v0','v1','vw'};

states1 = {'U','Lf','Ls','Lfp','Lsp','Isc','I','E','Rlo','Rhi','R','Dx','Tx'};

[i, s, d, lim] = get_addresses({states1, gps.vacc, gps.hiv}, [], [], [], 0);
d = char(d);

s.infectious  = [s.Isc, s.I, s.E, s.Dx];
s.symptomatic = [s.I, s.E, s.Dx];
s.prevalent   = [s.infectious, s.Tx];


% --- Include the auxiliaries ---------------------------------------------
names = {'inc','noti','mort','Dx'};
lgths =     [3,     3,     2,   1];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    i.aux.(names{ii}) = inds;
    lim = inds(end);
end
i.nx = lim;


% --- Make aggregators and selectors --------------------------------------

% Selectors for the incidence
tmp = zeros(3,i.nstates);
tmp(1,intersect(s.Isc,s.h0)) = 1;
tmp(2,intersect(s.Isc,s.h1)) = 1;
tmp(3,intersect(s.Isc,s.hart)) = 1;
agg.inc = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.Isc,:) = 1;
tmp(s.h1,s.h0) = 0;   tmp(s.h0,s.h1) = 0; 
tmp(s.hart,s.h1) = 0; tmp(s.h1,s.hart) = 0;
tmp(s.v1,s.v0) = 0; tmp(s.vw,s.v1) = 0;
sel.inc = tmp - diag(diag(tmp));

% Selectors for notifications
tmp = zeros(3,i.nstates);
tmp(1,intersect(s.Tx,s.h0)) = 1;
tmp(2,intersect(s.Tx,s.h1)) = 1;
tmp(3,intersect(s.Tx,s.hart)) = 1;
agg.noti = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.Tx,:) = 1;
tmp(s.h1,s.h0) = 0;   tmp(s.h0,s.h1) = 0; 
tmp(s.hart,s.h1) = 0; tmp(s.h1,s.hart) = 0;
tmp(s.v1,s.v0) = 0; tmp(s.vw,s.v1) = 0;
sel.noti = tmp - diag(diag(tmp));

% Selectors for counting presentations to care
tmp = zeros(1,i.nstates);
tmp(1,s.Dx) = 1;
agg.Dx = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.Dx,:) = 1;
tmp(s.h1,s.h0) = 0;   tmp(s.h0,s.h1) = 0; 
tmp(s.hart,s.h1) = 0; tmp(s.h1,s.hart) = 0;
tmp(s.v1,s.v0) = 0; tmp(s.vw,s.v1) = 0;
sel.Dx = tmp - diag(diag(tmp));
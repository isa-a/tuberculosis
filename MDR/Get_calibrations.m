clear all;

gps.strains = {'DS','MDR'};
gps.sectors = {'pu','pr'};

states0 = {'U'};                                                           % Unstructured
states1 = {'L','I','E','Rlo','Rhi','R'};                                   % + structured by strain
states2 = {'Dx','Tx','Tx2'};                                               % + structured by sector

[i, s, d, lim] = get_addresses({states0}, [], [], [], 0);
[i, s, d, lim] = get_addresses({states1, gps.strains}, i, s, d, lim);
[i, s, d, lim] = get_addresses({states2, gps.strains, gps.sectors}, i, s, d, lim);
d = char(d);
% Include the auxiliaries
i.aux.inc = i.nstates + [1:3];
i.nx = i.aux.inc(end);

s.infectious = [s.I, s.Dx, s.E, intersect(s.Tx, s.MDR)];
s.prevalent  = unique([s.infectious, s.Tx, s.Tx2]);

return;

% Selectors for the incidence
tmp = zeros(2,i.nstates); 
tmp(1,intersect(s.I,s.DS)) = 1;
tmp(2,intersect(s.I,s.MDR)) = 1;
agg.inc = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.I,:) = 1;
sel.inc = tmp - diag(diag(tmp));

tmp = zeros(i.nstates);
tmp(intersect(s.Tx,s.MDR),intersect(s.Tx,s.DS)) = 1;
sel.acqu = tmp - diag(diag(tmp));

% Natural history
r.reactivation  = 0.001;
r.careseeking2  = 12;
r.self_cure = 1/6;
r.relapse   = [0.032 0.14 0.0015];                                         % Recurrence rates for: good treatment; defaulters; non-adherers; combined failures
r.mort      = 1/66;
r.mort_TB   = 1/6;
p.Fast      = 0.14;
p.imm       = 0.8;

% Diagnosis and linkage to treatment
p.pu       = 0.5;
r.Dx       = 52;
p.Dx       = [0.83 0.75];
p.Tx_init  = [0.88 0.75];
r.MDR_acqu = 0.01;

p.MDR_rec  = [0.12 0];
p.Tx_init2 = [0.88 0];
p.SL_trans = [0.88 0];

% FL treatment
r.Tx      = 2;
pdef      = [.15 0.6];
r.default = r.Tx*pdef./(1-pdef);
p.cure    = [1 1];

% SL treatment
r.Tx2      = 0.5;
pdef       = [.15 0];
r.default2 = r.Tx2*pdef./(1-pdef);
p.cure2    = [0.5 0];

% Interventions
r.access  = 0;


p.kappa   = 2;                          % Relative infectiousness, post careseeking vs pre

prm.p = p; prm.r = r; 
ref.i = i; ref.s = s; ref.d = d;

data.prev = 160/.6;
data.incd = 217;
data.pmdr = 0.04;


obj = @(x) get_objective_wRNTCP(x, prm, ref, sel, agg, gps, data);
x1 = fminsearchbnd(obj, [6 6 1], [0 0 0], [], optimset('Display','iter'))
[out, aux] = get_objective_wRNTCP(x1, prm, ref, sel, agg, gps, data);

save('estims.mat');

% Show the recent incidence trajectories 
sol = aux.soln; t = sol(:,end);
inc = diff(interp1(t,sol(:,i.aux.inc),t(1):t(end)),1);
figure; 
subplot(1,2,1); plot(t(1):t(end-1),sum(inc(:,2:3)*1e5,2));
yl = ylim; yl(1) = 0; ylim(yl);
subplot(1,2,2); plot(t(1):t(end-1),inc(:,1)*1e5);
yl = ylim; yl(1) = 0; ylim(yl);
% What is the mean annual rate of TB decline
sel = sum(inc(:,[1,2]),2); dec_rate = -log(sel(end)/sel(1))/(t(end-1)-t(1))*100
% wRNTCP: version to include scale-up of RNTCP services

function [out, aux] = get_objective_wRNTCP(x, prm, ref, sel, agg, gps, data)

r = prm.r; p = prm.p; i = ref.i; s = ref.s;

r.beta = x(1:2);
r.careseeking = x(3);

% --- Set up the necessary models -----------------------------------------

% Final conditions
M2 = make_model(p, r, i, s, gps);

% In absence of RNTCP
p1 = p; r1 = r;
p1.pu = 0;
M1 = make_model(p1, r1, i, s, gps);

% In absence of MDR and RNTCP
p0 = p1; r0 = r1;
p0.pu = 0; r0.beta(2) = 0; r0.MDR_acqu = 0;
M0 = make_model(p0, r0, i, s, gps);


% --- Solve the models ----------------------------------------------------

% Equilibrium model
init = zeros(1,i.nx); seed = 1e-6; init(i.U) = (1-seed); init(i.I.DS) = seed;
geq = @(t,in)goveqs_basis2(t, in, M0, i, s, p0, sel, agg);
[t0, soln0] = ode15s(geq, [0 2e3], init, odeset('NonNegative',[1:i.nstates]));

% Introduce MDR-TB
init = soln0(end,:);
geq = @(t,in)goveqs_basis2(t, in, M1, i, s, p1, sel, agg);
[t1, soln1] = ode15s(geq, [1970 1997], init, odeset('NonNegative',[1:i.nstates]));

% RNTCP scale-up
init = soln1(end,:);
[t2, soln2] = ode15s(@(t,in) goveqs_scaleup(t, in, M1, M2, [1997 2007], i, s, p, sel, agg), [1997 2016], init, odeset('NonNegative',[1:i.nstates]));

soln = [soln2, t2];
sfin = soln(end,:);


% --- Get the objectives --------------------------------------------------

prev = sum(sfin(s.prevalent))*1e5;
tmp  = diff(interp1(t2,soln2(:,i.aux.inc),t2(1):t2(end)),1);
tmp = tmp(end,:)*1e5; incd = sum(tmp(1:2));
% pmdr = tmp(2)/incd;
pmdr = sum(tmp(2:3))/incd;

% Compose the objective function
sqdiff = @(dat, sim) sum((1-sim./dat).^2);
out = sqdiff(prev, data.prev) + sqdiff(incd, data.incd) + sqdiff(pmdr, data.pmdr);
out = out*100;

% --- Get additional outputs and package ----------------------------------

aux.soln  = [soln2, t2];
aux.soln0 = [soln0, t0];
aux.prev  = prev;
aux.incd  = incd;
aux.pmdr  = pmdr;
aux.M2    = M2;

end

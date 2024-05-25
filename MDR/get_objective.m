function [out, aux] = get_objective(x, prm, ref, gps, data)

r = prm.r; p = prm.p;
i = ref.i; s = ref.s;

r.beta = x(1:2);
r.careseeking = x(3:4);

M = make_model(p, r, i, s, gps);

init = zeros(1,i.nstates); seed = 1e-6; 
init(i.U.urb) = p.urb*(1-seed); init(i.I.urb) = p.urb*seed; 
init(i.U.rur) = (1-p.urb)*(1-seed); init(i.I.rur) = (1-p.urb)*seed; 
geq = @(t,in)goveqs_basis2(t, in, M, i, s, p);
[t, soln] = ode15s(geq, [0 1e3], init, odeset('NonNegative',[1:i.nstates]));

sfin = soln(end,:);

purb = sum(sfin(s.urb)); prur = sum(sfin(s.rur));
mat = M.lambda; mat(:,s.urb) = mat(:,s.urb)/purb; mat(:,s.rur) = mat(:,s.rur)/prur; 
arti = mat*sfin'*100;
% [tmp1, tmp2] = goveqs_basis(0, sfin', M, i, p); arti = tmp2*100;


prev = [sum(sfin(intersect(s.prevalent,s.urb))), sum(sfin(intersect(s.prevalent,s.rur)))]*1e5;
prev = prev./[purb, prur];

% Compose the objective function
sqdiff = @(dat, sim) sum((1-sim./dat).^2);
out = sqdiff(arti', data.arti) + sqdiff(prev, data.prev);
out = out*100;

aux.soln = [soln, t];
aux.arti = arti;
aux.prev = prev;
aux.M = M;

% save tempor;
% save temp;
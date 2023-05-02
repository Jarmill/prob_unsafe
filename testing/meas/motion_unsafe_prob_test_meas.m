%find the minimum distance between trajectories of the 'flow' system and a
%half-circle unsafe set

rng(343, 'twister');


%options 
SOLVE = 1;

n = 2;
order = 5;
d = 2*order; %0.2793, f = [0; 0], g=[0; 0.1];


%% problem parameters
% f_func = @(x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2) ];
% f_func = @(x) [0; 0]; %order 5: prob = 0.2793
f_func = @(x) [-x(2); x(1)]; %order 5: prob = 0.0877
% g_func = @(x) [1; 1]*0.4;
g_func = @(x) [0; 1]*0.1;

% Tmax = 5;
Tmax = 2;
BOX = 1.5;

%initial set
X0 = [0; 0.75];
% X0 = [0; 0.5];


%unsafe set
theta_c = 5*pi/4;  
% Cu = [-0.5; -0.75];
Cu = [0; 0];
Ru = 0.5;

%plotting

if SOLVE
mset clear
mset('yalmip',true);
mset(sdpsettings('solver', 'mosek'));
mpol('x', 2, 1);

f = Tmax*f_func(x);
g = Tmax*g_func(x);

% f = Tmax*f_func(x);
%initial set
%C0 = [1.2; 0];

% X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);

%half-circle set

% Cu = [0; -0.5];

c1f = Ru^2 - (x(1) - Cu(1)).^2 - (x(2) - Cu(2)).^2;

% theta_c = 3*pi/2;
% theta_c = 
w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(x(1) - Cu(1)) + w_c(2) * (x(2) - Cu(2)); 




%% set up measures
mpol('t0', 1, 1);
mpol('x0', 2, 1);
mu0 = meas([t0; x0]);

mpol('t_occ', 1, 1);
mpol('x_occ', 2, 1);
mu_occ = meas([t_occ; x_occ]);

mpol('tp', 1, 1);
mpol('xp', 2, 1);
mup = meas([tp; xp]);

mpol('tc', 1, 1);
mpol('xc', 2, 1);
muc = meas([tc; xc]);
%wasserstein

% u_cons = subs_vars([c1f; c2f], x, xp);
u_cons = subs_vars([c1f], x, xp);

% X0_con = subs_vars((x(1)-C0(1))^2 + (x(2)-C0(2))^2, x, x0) <= R0^2;
% X0_con = (x0 == C0);

%% support constraints
supp_con = [t0 == 0; x0 == X0;
    t_occ*(1-t_occ)>=0; tp*(1-tp) >=0; tc*(1-tc)>=0
    x_occ.^2 <= BOX^2; xp.^2 <= BOX^2; xc.^2 <= BOX^2;
    u_cons >= 0;
    ];



%% moment constraints

%get to writing

%liouville constraint
y0 = mom(mmon([t0; x0], d));
yp = mom(mmon([tp; xp], d));
yc = mom(mmon([tc; xc], d));
MOM_SUBS =0;

v  = mmon([t_occ; x_occ], d);
f_occ = subs_vars(f, x, x_occ);
g_occ = subs_vars(g, x, x_occ);
Ay_f = mom(diff(v, t_occ) + diff(v, x_occ)*f_occ); 
Ay_hess = Ay_f;
for i = 1:length(Ay_hess)
    h_curr = diff(diff(v(i), x_occ)', x_occ);
    Ay_hess(i) = mom(0.5*g_occ'*h_curr*g_occ);
end

Ay = Ay_f+ Ay_hess;
if MOM_SUBS
    Liou_con = yc == y0 + Ay - yp;
else
    Liou = Ay + (y0 - yp - yc);
    Liou_con = Liou == 0;
end

mom_con = [Liou_con; mass(mu0)==1];


%objective
objective = max(mass(mup));
% objective = min(mass(mu0));

%% solve problem
%Input LMI moment problem
P = msdp(objective, ...
    mom_con, supp_con);

%solve LMIP moment problem
[status, obj, m, dual_rec] = msol(P);    

%% analyze solutions
prob_unsafe = double(mass(mup));

fprintf('Prob Unsafe<= %0.4f\n', prob_unsafe)
M0 = double(mmat(mu0));
Mocc = double(mmat(mu_occ));
Mp = double(mmat(mup));
Mc = double(mmat(muc));

M0_1 = M0(1:(n+2), 1:(n+2));
Mp_1 = Mp(1:(n+2), 1:(n+2));
Mc_1 = Mc(1:(2*n+1), 1:(2*n+1));

rankp = rank(Mp_1, 1e-3);
rank0 = rank(M0_1, 1e-3);
rankw = rank(Mc_1, 1e-3);

% xu_rec = double(mom(xu));
xp_rec = double(mom(xp));
x0_rec = double(mom(x0));
tp_rec = Tmax*double(mom(tp));


% optimal_pt = all([rankp; rank0; rankw]==1);

end 

%floating platform dynamics https://arxiv.org/pdf/2208.10752
%% variables and dynamics
t = sdpvar(1,1);
x = sdpvar(4,1);

PLOT = 1;

% b = -0.1;
sigma = 0.05;
f = f_dyn(t, x);
g =  sigma * [0; -x(2)+1; 0; -x(4)+1];


%% unsafe set

%% put together the constraints

x0 = [-0.25; 0.2; 0.5; 0.1];


%% Support Sets
T = 8;
% Xmax = 1.5;
Xmax = 1.5;
Vmax = 1;
Ru = 0.15;
% Xmax = 1.25;
% Xmax = 2.5;
Xbox = [Xmax.^2-x([1, 3]).^2; Vmax.^2-x(1+[1, 3]).^2];
X = struct('ineq', Xbox, 'eq', []);
Xall = struct('ineq', [t*(1-t); Xbox], 'eq', []);

unsafe_cons_1 = Ru^2 - (x(3)-x(1))^2;
Xu1 = struct('ineq', unsafe_cons_1, 'eq', []);
Xu1all = struct('ineq', [t*(1-t); Xu1.ineq; Xbox], 'eq', []);

%% polynomials
%linear system
% % order = 1; %1.000
% % order = 2; %0.4371
% % order = 3; %0.1732
% % order = 4; %0.1137
% order = 5;   %0.1118

%nonlinear spring
% order = 1; %1.000
% order = 2; %   0.9236
% order = 3; %    0.3068
order = 4; %
% order = 5;   %

d = 2*order;

[v, cv, mv] = polynomial([t;x], d);
[w, cw, mw] = polynomial([x], d);

fT = T*f;
gT = T*g;

Lv = jacobian(v, t) + jacobian(v, x)*fT + 0.5*(gT)'*hessian(v, x)*(gT);

v0 = replace(v, [t], [0]);

%toggles for different experiments
LEBESGUE = 0;
CIRC = 0;

consinit = [];
coeffinit = [];

% if LEBESGUE
%     leb = LebesgueBoxMom( d, [-Xmax, Xmax; -Xmax, Xmax]',1);
%     w_leb = cw'*leb;
%     objective = w_leb;
% else
    if CIRC
        gamma = sdpvar(1, 1);
        objective = gamma;
        R0 = 0.2;
        X0 = struct('ineq', R0^2 - sum((x-x0).^2), 'eq', []);
        [put_init, consinit, coeffinit] = constraint_psatz(gamma-v0, X0, x, d);

    else
        v0x = replace(v0, x, x0);
        objective = v0x;
    end
% end

zero_con = (coefficients(v0 - w, x)==0);

[put_lie, conslie, coefflie] = constraint_psatz(-Lv, Xall, [t;x], d);

[put_unsafe1, consunsafe1, coeffunsafe1] = constraint_psatz(v-1, Xu1all, [t;x], d);

consunsafe = [consunsafe1];
coeffunsafe = [coeffunsafe1];

%% exit measure
%this is the complement measure
%decompose it into parts
vT = replace(v, t, 1);

[put_base_T, consbase, coeffbase] = constraint_psatz(v, Xall, [t;x], d);


%% assemble all constraints

cons = [conslie:'lie'; consbase:'base'; consunsafe:'unsafe'; zero_con:'zero'; consinit];
coeff = [coefflie; coeffbase; coeffunsafe; cv; cw; coeffinit];

opts = sdpsettings('solver', 'mosek');


[sol,u,Q] = solvesos(cons,objective,opts,coeff);

disp(value(objective))
% recovery
v_rec = value(cv)'*mv;
v0_rec = replace(v_rec, t, 0);
obj_rec = value(objective);

%% plotting
% 
% if PLOT
% wp = polyval_func(v0_rec, x);
% figure(2)
% clf

%% visualize

function dx = f_dyn(t, x)
%dynamics for transient stability
x1 = x(1:2);
x2 = x(3:4);

k12 = 1;
k12_c = -0.05;
d12 = 0.5;

xdiff = ((x2(1)-x1(1))-d12);

f12 = k12*xdiff + k12_c*xdiff^3;
% f12 = k*(-x1(1)-d);
% f23 = k*(x3(1)-d);

vd1 = f12 - x1(2);
vd2 = -x2(1) - x2(2) -f12 ;

dx = [x1(2); vd1; x2(2); vd2];

end
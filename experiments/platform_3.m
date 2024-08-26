%2d motion when drift term is cons\tant

%% variables and dynamics
t = sdpvar(1,1);
x = sdpvar(6,1);

PLOT = 1;

% b = -0.1;
sigma = 0.1;
f = f_dyn(t, x);
g =  sigma * [0; -x(2)+1; 0; -x(4)+1; 0; -x(6)+1];


%% unsafe set

%% put together the constraints

x0 = [-1; 0.3; 0.5; 0.3; 1; -0.1];


%% Support Sets
T = 8;
% Xmax = 1.5;
Xmax = 2;
Vmax = 1;
% Xmax = 1.25;
% Xmax = 2.5;
Xbox = [Xmax.^2-x([1, 3, 5]).^2; Vmax.^2-x(1+[1, 3, 5]).^2];
X = struct('ineq', Xbox, 'eq', []);
Xall = struct('ineq', [t*(1-t); Xbox], 'eq', []);


unsafe_cons_1 = (x(3)-x(1))^2 - 0.2;
unsafe_cons_2 = (x(5)-x(3))^2 - 0.2;
Xu1 = struct('ineq', unsafe_cons_1, 'eq', []);
Xu2 = struct('ineq', unsafe_cons_2, 'eq', []);

Xu1all = struct('ineq', [t*(1-t); Xu1.ineq; Xbox], 'eq', []);
Xu2all = struct('ineq', [t*(1-t); Xu2.ineq; Xbox], 'eq', []);


%% polynomials
%polynomial definition
% order =6; 
order = 3;
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
[put_unsafe2, consunsafe2, coeffunsafe2] = constraint_psatz(v-1, Xu2all, [t;x], d);

consunsafe = [consunsafe1; consunsafe2];
coeffunsafe = [coeffunsafe1; coeffunsafe2];

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
x3 = x(5:6);

k12 = 1;
k23 = 0.75;
d12 = 0.5;
d23 = 0.7;

f12 = k12*((x2(1)-x1(1))-d12);
f23 = k23*((x3(1)-x2(1))-d23);

% f12 = k*(-x1(1)-d);
% f23 = k*(x3(1)-d);

vd1 = f12 - x1(2);
vd2 = -x2(1) - x2(2) -f12 + f23;
vd3 = -f23 - x3(2);

dx = [x1(2); vd1; x2(2); vd2; x3(2); vd3];

end
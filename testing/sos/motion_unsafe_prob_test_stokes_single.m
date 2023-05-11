%2d motion when drift term is constant

%% variables and dynamics
t = sdpvar(1,1);
x = sdpvar(2,1);

% f =  2*[-x(2); x(1)];
f = 2*[0; -1];
% g = 2*[0; 0.1];
g = [0; 0];



%% unsafe set

% theta_c = 5*pi/4; 
% Cu = [-0.5; -0.75]; %original parameters
% % Ru = 0.5;

Cu = [0; 0];
Ru = 0.5;

y = x;

c1f = Ru^2 - (y(1) - Cu(1)).^2 - (y(2) - Cu(2)).^2;

% unsafe_cons = [c1f; c2f];
unsafe_cons = c1f;

%% put together the constraints

% order = 1;
order = 2;
% order = 3;
% order = 4; 
% order = 5; 
% order = 6; 
d = 2*order;


x0 = [0; 0.75];

%% Support Sets
T = 2;
Xmax = 1.5;
Xall = struct('ineq', [t*(1-t); Xmax.^2-x.^2], 'eq', []);


Xu = struct('ineq', unsafe_cons, 'eq', []);
Xuall = struct('ineq', [t*(1-t); Xu.ineq; Xmax.^2-x.^2], 'eq', []);

%boundary
bXu = struct('ineq', [], 'eq', unsafe_cons);
% bXuall = struct('ineq', [t*(1-t)], 'eq', unsafe_cons);


unsafe_boundary = t*(1-t)*prod(unsafe_cons);
bXusingle= struct('ineq', [t*(1-t); unsafe_cons], 'eq', unsafe_boundary);

%% polynomials
%polynomial definition
[v, cv, mv] = polynomial([t;x], d);
[ut, cut, mut] = polynomial([t;x], d);
[ux1, cux1, mux1] = polynomial([t;x], d);
[ux2, cux2, mux2] = polynomial([t;x], d);

ux = [ux1; ux2];
cux = [cux1; cux2];
mux = [mux1; mux2];

fT = T*f;
gT = T*g;

Lv = jacobian(v, t) + jacobian(v, x)*fT + 0.5*(gT)'*hessian(v, x)*(gT);

% du = jacobian(ux1, x(1)) + jacobian(ux2, x(2)) + jacobian(ut, t);
ub = [ut; ux]*unsafe_boundary;
dub = jacobian(ub, [t; x]);
divub = trace(dub);

ut0 = replace(ut, t, 0);
ut1 = replace(ut, t, 1);

v0 = replace(v, [t;x], [0; x0]);

objective = v0;

dchange = 1;

[put_lie, conslie, coefflie] = constraint_psatz(-Lv, Xall, [t;x], d);
[put_base, consbase, coeffbase] = constraint_psatz(v, Xall, [t;x], d);
[put_unsafe, consunsafe, coeffunsafe] = constraint_psatz(v-1-divub, Xuall, [t;x], d+dchange);


cons = [conslie:'lie'; consbase:'base'; consunsafe:'unsafe'];
coeff = [coefflie; coeffbase; coeffunsafe; cv; cux; cut];

cons_nullu = [cux==0; cut==0];

NULL_U = 0;

if NULL_U   
    cons = [cons_nullu; cons];
end
opts = sdpsettings('solver', 'mosek');


[sol,u,Q] = solvesos(cons,objective,opts,coeff);


%% recovery
v_rec = value(cv)'*mv;
v0_rec = replace(v_rec, [t; x], [0; x0]);
obj_rec = value(objective);

fprintf('Prob unsafe <= %0.4f \n', obj_rec);
% disp(obj_rec)
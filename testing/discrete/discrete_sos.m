%2d discrete-time dynamics

%% variables and dynamics
t = sdpvar(1,1);
x = sdpvar(2,1);
omega = sdpvar(1, 1);

A = [-0.3, 0.8; -0.9, -0.1];
w0 = 0.5;
f = A*x - [0; 0.2*x(1)^2] + [(w0/2)*omega*prod(x); 0]; 

%% unsafe set

theta_c = pi/4; 
Cu = [0.8; 0.2];

% Cu = [0; 0];
Ru = 0.4;
y = x;

c1f = Ru^2 - (y(1) - Cu(1)).^2 - (y(2) - Cu(2)).^2;

% theta_c = 
w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(y(1) - Cu(1)) + w_c(2) * (y(2) - Cu(2)); 
% c2f = [];
unsafe_cons = [c1f; c2f];
% unsafe_cons = c1f;

%% put together the constraints

order = 2; 
d = 2*order;


% x0 = [0; 0.75];
x0=[-1; 0];

%% Support Sets
T = 10;
Xmax = 1.5;
Xall = struct('ineq', [t*(1-t); Xmax.^2-x.^2], 'eq', []);


Xu = struct('ineq', unsafe_cons, 'eq', []);
Xuall = struct('ineq', [t*(1-t); Xu.ineq; Xmax.^2-x.^2], 'eq', []);


%% polynomials
%polynomial definition
[v, cv, mv] = polynomial([t;x], d);
[w, cw, mw] = polynomial([x], d);


v_next = replace(v, [t; x], [t+1/T; f]);
[comega, momega] = coefficients(v_next, omega);
mom_omega_norm = NormalMom(d);
v_next_marg = comega'*mom_omega_norm;

%then marginalize out w;
Lv = v_next_marg - v;

v0 = replace(v, [t], [0]);

% objective = v0;

%toggles for different experiments
LEBESGUE = 0;
CIRC = 1;

consinit = [];
coeffinit = [];


if LEBESGUE
    leb = LebesgueBoxMom( d, [-Xmax, Xmax; -Xmax, Xmax]',1);
    w_leb = cw'*leb;
    objective = w_leb;
else
    if CIRC
        gamma = sdpvar(1, 1);
        objective = gamma;
        R0 = 0.4;
        X0 = struct('ineq', R0^2 - sum((x-x0).^2), 'eq', []);
        [put_init, consinit, coeffinit] = constraint_psatz(gamma-v0, X0, x, d);

    else
        v0x = replace(v0, x, x0);
        objective = v0x;
    end
end

zero_con = (coefficients(v0 - w, x)==0);


[put_lie, conslie, coefflie] = constraint_psatz(-Lv, Xall, [t;x], d);
[put_base, consbase, coeffbase] = constraint_psatz(v, Xall, [t;x], d);
[put_unsafe, consunsafe, coeffunsafe] = constraint_psatz(v-1, Xuall, [t;x], d);

cons = [conslie:'lie'; consbase:'base'; consunsafe:'unsafe'; zero_con:'zero'; consinit];
coeff = [coefflie; coeffbase; coeffunsafe; cv; cw; coeffinit];

opts = sdpsettings('solver', 'mosek');


[sol,u,Q] = solvesos(cons,objective,opts,coeff);


%% recovery
v_rec = value(cv)'*mv;
v0_rec = replace(v_rec, t, 0);
obj_rec = value(objective);

%% plotting
wp = polyval_func(v0_rec, x);
figure(1)
clf
fsurf(@(x, y) wp([x; y]), [-1.5, 1.5, -1.5, 1.5])
xlabel('x_1')
ylabel('x_2')
zlabel('Prob')
title(sprintf('Probability of Unsafety (T=%d)', T), 'FontSize', 16)


theta_half_range = linspace(theta_c-pi/2, theta_c + pi/2, 200);
circ_half = [cos(theta_half_range); sin(theta_half_range)];
Xu = Cu + circ_half* Ru;
view(-8, 65)

figure(2)
clf
hold on
fcontour(@(x, y) wp([x; y]), [-1.5, 1.5, -1.5, 1.5], 'levellist', 0:0.1:1, 'fill', 1)
patch(Xu(1, :), Xu(2, :), zeros(size(Xu(1, :))), 'r', 'Linewidth', 3, 'EdgeColor', 'none')
axis square
xlabel('x_1')
ylabel('x_2')
colorbar
title(sprintf('Probability of Unsafety (T=%d)', T), 'FontSize', 16)
% fprintf('Prob unsafe <= %0.4f \n', obj_rec);
% disp(obj_rec)

%2d motion when drift term is constant

%% variables and dynamics
t = sdpvar(1,1);
x = sdpvar(2,1);

% b = -0.1;
sigma = 0.1;
% f =  [x(2); -(x(1) +x(2) + 0.5*x(1)^3)];
g = sigma * [0;1];
f =  [-x(2); x(1)];
% f = [0; 0];
% f = 2*[0; -1];
% g = 2*[0; 0.1];

%% unsafe set

theta_c = 5*pi/4; 
Cu = [-0.5; -0.75]; %original parameters
% % Ru = 0.5;

% Cu = [0; 0];
Ru = 0.5;

y = x;

c1f = Ru^2 - (y(1) - Cu(1)).^2 - (y(2) - Cu(2)).^2;

theta_c = 3*pi/2;
% theta_c = 
w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(y(1) - Cu(1)) + w_c(2) * (y(2) - Cu(2)); 
% c2f = [];
unsafe_cons = [c1f; c2f];
% unsafe_cons = c1f;

%% put together the constraints

order = 6; 
d = 2*order;


x0 = [0; 0.75];

%% Support Sets
T = 1;
% T = 3;
% T = 5;
% T = 8;
Xmax = 1.5;
Xall = struct('ineq', [t*(1-t); Xmax.^2-x.^2], 'eq', []);


Xu = struct('ineq', unsafe_cons, 'eq', []);
Xuall = struct('ineq', [t*(1-t); Xu.ineq; Xmax.^2-x.^2], 'eq', []);


%% polynomials
%polynomial definition
[v, cv, mv] = polynomial([t;x], d);
[w, cw, mw] = polynomial([x], d);

fT = T*f;
gT = T*g;

Lv = jacobian(v, t) + jacobian(v, x)*fT + 0.5*(gT)'*hessian(v, x)*(gT);

v0 = replace(v, [t], [0]);

% objective = v0;

leb = LebesgueBoxMom( d, [-Xmax, Xmax; -Xmax, Xmax]',1);
w_leb = cw'*leb;
objective = w_leb;

zero_con = (coefficients(v0 - w, x)==0);

[put_lie, conslie, coefflie] = constraint_psatz(-Lv, Xall, [t;x], d);
[put_base, consbase, coeffbase] = constraint_psatz(v, Xall, [t;x], d);
[put_unsafe, consunsafe, coeffunsafe] = constraint_psatz(v-1, Xuall, [t;x], d);

cons = [conslie:'lie'; consbase:'base'; consunsafe:'unsafe'; zero_con:'zero'];
coeff = [coefflie; coeffbase; coeffunsafe; cv; cw];

opts = sdpsettings('solver', 'mosek');


[sol,u,Q] = solvesos(cons,objective,opts,coeff);


%% recovery
v_rec = value(cv)'*mv;
v0_rec = replace(v_rec, t, 0);
obj_rec = value(objective);

%% plotting
wp = polyval_func(v0_rec, x);
figure(1)
fsurf(@(x, y) wp([x; y]), [-1.5, 1.5, -1.5, 1.5])
xlabel('x_1')
ylabel('x_2')
zlabel('Prob')
title(sprintf('Unsafety Upper-Bound (T=%d)', T), 'FontSize', 16)


theta_half_range = linspace(theta_c-pi/2, theta_c + pi/2, 200);
circ_half = [cos(theta_half_range); sin(theta_half_range)];
Xu = Cu + circ_half* Ru;
view(-8, 65)

figure(2)
hold on
fcontour(@(x, y) wp([x; y]), [-1.5, 1.5, -1.5, 1.5], 'levellist', 0:0.1:1, 'fill', 1)
patch(Xu(1, :), Xu(2, :), zeros(size(Xu(1, :))), 'r', 'Linewidth', 3, 'EdgeColor', 'none')
axis square
xlabel('x_1')
ylabel('x_2')
colorbar
title(sprintf('Unsafety Upper-Bound (T=%d)', T), 'FontSize', 16)
% fprintf('Prob unsafe <= %0.4f \n', obj_rec);
% disp(obj_rec)

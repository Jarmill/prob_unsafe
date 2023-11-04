%2d motion when drift term is constant

%% variables and dynamics
t = sdpvar(1,1);
x = sdpvar(2,1);

PLOT = 1;

% b = -0.1;
sigma = 0.1;
f =  [x(2); -(x(1) +x(2) + 0.5*x(1)^3)];
g = sigma * [0;1];
% f =  [-x(2); x(1)];
% f = [0; 0];
% f = 2*[0; -1];
% g = 2*[0; 0.1];

%% unsafe set

%% unsafe set (moon)
%plot the moon
h_in = 0.3;
h_out = 1;
 
y = x;
% moon_center = [0.4;0.0];
% moon_theta = -3*pi/2;
% moon_scale = 0.5;

moon_center = [0.2;0.0];
moon_theta = -11*pi/8;
moon_scale = 0.5;


%statistics of the moon
c_in = [0;0.5*(1/h_in - h_in)];
r_in = 0.5*(1/h_in + h_in);

c_out = [0;0.5*(1/h_out - h_out)];
r_out = 0.5*(1/h_out + h_out);

c_in_scale = moon_rot*c_in*moon_scale + moon_center;
c_out_scale = moon_rot*c_out*moon_scale + moon_center;

r_in_scale = moon_scale*r_in;
r_out_scale = moon_scale*r_out;

%constraints of the moon
con_inner =  sum((y-c_in_scale).^2) - r_in_scale^2;
con_outer =  -sum((y-c_out_scale).^2) + r_out_scale^2;

unsafe_cons = [con_inner; con_outer];

% unsafe_cons = c1f;

%% put together the constraints

% x0 = [0; 0.75];
% x0 = [0.75; 0];
x0 = [0.85; -0.75]; %order 6 prb <= 3.0565354310e-01 (need to simulate)



%% Support Sets
T = 5;
% Xmax = 1.5;
% Xmax = 2;
Xmax = 1.25;
% Xmax = 2.5;
Xbox = Xmax.^2-x.^2;
Xall = struct('ineq', [t*(1-t); Xbox], 'eq', []);


Xu = struct('ineq', unsafe_cons, 'eq', []);
Xuall = struct('ineq', [t*(1-t); Xu.ineq; Xbox], 'eq', []);


%% polynomials
%polynomial definition
order =6; 
% order = 4;
d = 2*order;

[v, cv, mv] = polynomial([t;x], d);
[w, cw, mw] = polynomial([x], d);

fT = T*f;
gT = T*g;

Lv = jacobian(v, t) + jacobian(v, x)*fT + 0.5*(gT)'*hessian(v, x)*(gT);

v0 = replace(v, [t], [0]);

%toggles for different experiments
LEBESGUE = 1;
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
        R0 = 0.2;
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

disp(value(objective))
% recovery
v_rec = value(cv)'*mv;
v0_rec = replace(v_rec, t, 0);
obj_rec = value(objective);

%% plotting

if PLOT
wp = polyval_func(v0_rec, x);
figure(2)
clf

fsurf(@(x, y) min(1, wp([x; y])), [-1.5, 1.5, -1.5, 1.5])
xlabel('x_1')
ylabel('x_2')
xlim([-Xmax, Xmax])
ylim([-Xmax, Xmax])
zlabel('Prob')
zlim([0, 1])
title(sprintf('Unsafety Upper-Bound (T=%d)', T), 'FontSize', 16)


% theta_half_range = linspace(theta_c-pi/2, theta_c + pi/2, 200);
% circ_half = [cos(theta_half_range); sin(theta_half_range)];

%plot the moon
x_moon = moon_base(h_in, h_out);
moon_rot = [cos(moon_theta), sin(-moon_theta); sin(moon_theta), cos(moon_theta)];
x_moon_move = moon_rot*x_moon*moon_scale + moon_center;

Xu = x_moon_move;
% Xu = Cu + circ_half* Ru;
view(-8, 65)

%% visualize
figure(4)
clf
hold on
cont_tol = 0;
ll=[0:0.1:1, (1+cont_tol)];
fcontour(@(x, y) min(1+cont_tol, wp([x; y])), [-Xmax, Xmax, -Xmax, Xmax], 'levellist', ll, 'fill', 1, 'MeshDensity', 150)
patch(Xu(1, :), Xu(2, :), zeros(size(Xu(1, :))), 'r', 'Linewidth', 3, 'EdgeColor', 'none')

if ~LEBESGUE
    if CIRC
        theta_full = linspace(0, 2*pi, 200);        
%         circ = [cos(theta_full); sin(theta_full)];
        plot(R0*cos(theta_full) + x0(1), R0*sin(theta_full)+ x0(2), 'm', 'LineWidth', 3)
    else
        scatter(x0(1), x0(2), 200, 'magenta', 'filled')
    end
end
axis square
xlabel('x_1')
ylabel('x_2')
xlim([-Xmax, Xmax])
ylim([-Xmax, Xmax])
colorbar
title(sprintf('Unsafety Upper-Bound (T=%d)', T), 'FontSize', 16)
% fprintf('Prob unsafe <= %0.4f \n', obj_rec);
% disp(obj_rec)
end
%sample SDE for the motion system
rng(1, 'twister')


%dynamics
sigma = 0.1;
f =  @(t, x) [x(2); -(x(1) +x(2) + 0.5*x(1)^3)];
g =  @(t, x) sigma * [0;1];

%set geometry
Xmax = 2;
Ru = 0.5;
Cu = [-0.5; -0.75]; %original parameters
theta_c = 3*pi/2;
R0 = 0.2;
x0 = [0.85; -0.75];
dt = 1e-3;

param = struct('Ru', Ru,  'Cu', Cu, 'theta_c', theta_c, 'Xmax', Xmax, 'dt', dt);

Tmax = 5;

%sampler
x0 = [0.85; -0.75];
opts_fw = sdeset('EventsFun', @(t, x) stoch_event_motion(t, x, param),...
    'SDEType', 'Ito',  'DiagonalNoise', 'yes', 'ConstGFUN', 'yes');

%% perform sampling
Ntraj = 1000;
x_traj = cell(Ntraj, 1);
for i = 1:Ntraj
% x0_curr = x0;
x0_curr = R0*sphere_sample(1, 2)' + x0;
[x_traj{i},W,TE,YE,WE,IE] = sde_euler(f,g,0:param.dt:Tmax,...
        x0_curr,opts_fw);
end
xlabel('x_1')
ylabel('x_2')
title('Sampled Trajectories', 'FontSize', 16)


%% perform visualization
figure(3)
clf
hold on
%generate patches
theta_half_range = linspace(theta_c-pi/2, theta_c + pi/2, 200);
circ_half = [cos(theta_half_range); sin(theta_half_range)];
Xu = Cu + circ_half* Ru;

theta_full = linspace(0, 2*pi, 200);        
patch(Xu(1, :), Xu(2, :), zeros(size(Xu(1, :))), 'r', 'Linewidth', 3, 'EdgeColor', 'none')

for i = 1:Ntraj
    x_curr = x_traj{i};
    plot(x_curr(:, 1), x_curr(:, 2), 'c');
end

plot(R0*cos(theta_full) + x0(1), R0*sin(theta_full)+ x0(2), 'm', 'LineWidth', 3)



function [event_eval, terminal, direction] = stoch_event_motion(tp, xp, param)
    %event function for @ode15 or other solver
    %Start control of the system (it is too far away from the origin)
    Npt = size(xp, 2);
    event_eval = zeros(1, Npt);
    for i = 1:Npt
        x= xp(:, i);
        t= tp(:, i);               
        
        w_c = [cos(param.theta_c); sin(param.theta_c)];
        
%         crit_bk = sum(xcurr.^2);
        c1f = param.Ru^2 - (x(1) - param.Cu(1)).^2 - (x(2) - param.Cu(2)).^2;
        c2f = w_c(1)*(x(1) - param.Cu(1)) + w_c(2) * (x(2) - param.Cu(2)); 
        
        cbox = param.Xmax - abs(x);
        event_eval(i) = 2*all([c1f; c2f; cbox])-1; 
    end

    %stop integrating when the system falls outside support

    terminal = 1;
    direction = 0;                        
end
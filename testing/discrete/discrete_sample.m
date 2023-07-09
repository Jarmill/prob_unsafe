rng(23, 'twister');

%dynamical system parameters

% A = [-0.3, 0.8; -0.75, -0.3];
A = [-0.3, 0.8; -0.9, -0.1];
% A = [-0.3, 0.8; -0.75, -0.];
% knn = 0.25;
% + [knn*x(1)*x(2); 0];
f = @(x) A*x - [0; 0.2*x(1)^2]; 
% C0=[-1.5; 0];
C0=[-1; 0];
R0 = 0.4;
% w0 = 0.2;
% w0 = 0.4;
w0 = 0.5;
% w0 = 1;


%% unsafe set

% theta_c = 5*pi/4; 
theta_c = pi/4; 
% Cu = [-0.8; -0.8]; %original parameters
Cu = [0.8; 0.2];
% % Ru = 0.5;

% Cu = [0; 0];
Ru = 0.4;


%% sample
T = 10;

smp = struct('x', @() R0*ball_sample(1, 2)' + C0, 'w', @() randn(1)*w0/2);
% smp = struct('x', @() R0*sphere_sample(1, 2)' + C0, 'w', @() (2*rand(1)-1) * w0);

x_curr = x0;
Ntraj = 1500;
traj = cell(Ntraj, 1);
for k=1:Ntraj
    traj{k} = curr_traj(smp, T, f);
end  

%% visualize

figure(1)
clf
hold on
theta_half_range = linspace(theta_c-pi/2, theta_c + pi/2, 200);
circ_half = [cos(theta_half_range); sin(theta_half_range)];
Xu = Cu + circ_half* Ru;
patch(Xu(1, :), Xu(2, :), zeros(size(Xu(1, :))), 'r', 'Linewidth', 3, 'EdgeColor', 'none')

for k=1:Ntraj
    ct = traj{k};
    plot(ct(1, :), ct(2, :), '.c')
end  

th = linspace(0, 2*pi, 300);
plot(R0*cos(th)+C0(1), R0*sin(th)+C0(2), 'k', 'linewidth', 2)
pbaspect([diff(xlim), diff(ylim), 1])


% axis square
xlabel('x_1')
ylabel('x_2')
title('Sampled Trajectories', 'Fontsize', 16)

%% function
function x= curr_traj(smp, T, f)

x0 = smp.x();
x = zeros(2, T+1);
x(:, 1) = x0;

x_curr = x0;

for i = 1:T
    wc = smp.w();
%     x_next = f(x_curr) + [0; wc];
    x_next = f(x_curr) + [wc*prod(x_curr); 0];
    x(:, i+1) = x_next;
    x_curr = x_next;
end
end
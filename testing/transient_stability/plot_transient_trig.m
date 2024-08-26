%sample SDE for transient stability system
%requires https://www.mathworks.com/matlabcentral/fileexchange/56406-sdetools
% https://arxiv.org/pdf/1811.01372
rng(1, 'twister')


%dynamics
sigma = 0.15;
f =  @(t, x) f_dyn(t, x);
g =  @(t, x) sigma * [1; 1; 0; 0; 0; 0];

%set geometry
dt = 1e-3;

% param = struct('Ru', Ru,  'Cu', Cu, 'theta_c', theta_c, 'Xmax', Xmax, 'dt', dt);
param = struct('dt', dt);
Tmax = 8;
% Tmax = 30;
R0 = 0.15;
Xmax = 0.8;


%sampler
% x0 = [0.02; 0.06; 0; 0];
x0 = [0; 0; 0; 0];
opts_fw = sdeset('EventsFun', @(t, x) stoch_event_motion(t, x, param),...
    'SDEType', 'Ito',  'DiagonalNoise', 'yes', 'ConstGFUN', 'yes');

synclim = 0.2;


%% perform sampling
Ntraj = 300;
x_traj = cell(Ntraj, 1);
for i = 1:Ntraj
% x0_curr = x0;
x0_curr = [R0*(2*rand(2, 1) - 1); 0; 0] + x0;
[x_traj{i},W,TE,YE,WE,IE] = sde_euler(f,g,0:param.dt:Tmax,...
        x0_curr,opts_fw);
end

%% plot
figure(3)
clf
hold on
%generate patches
% theta_half_range = linspace(theta_c-pi/2, theta_c + pi/2, 200);
% circ_half = [cos(theta_half_range); sin(theta_half_range)];
% Xu = Cu + circ_half* Ru;
% 
% theta_full = linspace(0, 2*pi, 200);        
% patch(Xu(1, :), Xu(2, :), zeros(size(Xu(1, :))), 'r', 'Linewidth', 3, 'EdgeColor', 'none')

for i = 1:Ntraj
    x_curr = x_traj{i};
    plot(x_curr(:, 1), x_curr(:, 2), 'c');
end
xlabel('\theta_1')
ylabel('\theta_2')
axis square
xl = [-Xmax, Xmax];
yl = [-Xmax, Xmax];

patch(R0*[-1, -1, 1, 1, -1], R0*[-1, 1, 1, -1, -1], 'm', 'EdgeColor', 'none');

fcontour(@(x, y) 1 - cos(x).*cos(y) - sin(x).*sin(y), [xl, yl], 'levellist', [synclim], 'linecolor', 'r')
xlim(xl);
ylim(yl);



% title('Sampled Trajectories', 'FontSize', 16)

figure(5)
clf
hold on
for i = 1:Ntraj
    x_curr = x_traj{i};
    % plot((cos(x_curr(:, 1))- cos(x_curr(:, 2))).^2 + (sin(x_curr(:, 1))- sin(x_curr(:, 2))).^2, 'c');

plot(1- (x_curr(:, 3).*x_curr(:, 4) + x_curr(:, 5).*x_curr(:, 6)), 'c');

end


function dx = f_dyn(t, x)
%dynamics for transient stability
w = x(1:2);
c = x(3:4);
s = x(5:6);


sin_diff = s(1)*c(2) - c(1)*s(2);
fw1 = -s(1) -0.4*w(1) -0.5*(sin_diff);
fw2 = -0.5*s[2] - 0.5*(-sin_diff) - 0.5*w[2] + 0.05;
fw = [fw1; fw2];
f = [fw; s.*w; (-c).*w];

dx = f;
end

function [event_eval, terminal, direction] = stoch_event_motion(tp, xp, param)
    %event function for @ode15 or other solver
    %Start control of the system (it is too far away from the origin)
%     Npt = size(xp, 2);
%     event_eval = zeros(1, Npt);
%     for i = 1:Npt
%         x= xp(:, i);
%         t= tp(:, i);               
% 
%         w_c = [cos(param.theta_c); sin(param.theta_c)];
% 
% %         crit_bk = sum(xcurr.^2);
%         c1f = param.Ru^2 - (x(1) - param.Cu(1)).^2 - (x(2) - param.Cu(2)).^2;
%         c2f = w_c(1)*(x(1) - param.Cu(1)) + w_c(2) * (x(2) - param.Cu(2)); 
% 
%         cbox = param.Xmax - abs(x);
%         event_eval(i) = 2*all([c1f; c2f; cbox])-1; 
%     end

    event_eval = 0;

    %stop integrating when the system falls outside support

    terminal = 1;
    direction = 0;                        
end
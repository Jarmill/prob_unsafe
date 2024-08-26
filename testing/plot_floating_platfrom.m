%sample SDE for transient stability system
%requires https://www.mathworks.com/matlabcentral/fileexchange/56406-sdetools
% https://arxiv.org/pdf/1811.01372
rng(1, 'twister')


%dynamics
% sigma = 0.25;
sigma = 0.05;
f =  @(t, x) f_dyn(t, x);
g =  @(t, x) sigma * [0; -x(2)+1; 0; -x(4)+1; 0; -x(6)+1];

%set geometry
dt = 1e-3;

% param = struct('Ru', Ru,  'Cu', Cu, 'theta_c', theta_c, 'Xmax', Xmax, 'dt', dt);
param = struct('dt', dt);
Tmax = 8;
% Tmax = 30;


%sampler
% x0 = [-1; 0; 0; 0; 1; 0];
x0 = [-1; 0.3; 0.5; 0.3; 1; -0.1];
opts_fw = sdeset('EventsFun', @(t, x) stoch_event_motion(t, x, param),...
    'SDEType', 'Ito',  'DiagonalNoise', 'yes', 'ConstGFUN', 'yes');



%% perform sampling
Ntraj = 100;
x_traj = cell(Ntraj, 1);
for i = 1:Ntraj
x0_curr = x0;
% x0_curr = [R0*(2*rand(2, 1) - 1); 0; 0] + x0;
[x_traj{i},W,TE,YE,WE,IE] = sde_euler(f,g,0:param.dt:Tmax,...
        x0_curr,opts_fw);
end

%% plot
figure(3)
clf
hold on

for i = 1:Ntraj
    x_curr = x_traj{i};
    plot(x_curr(:, 1), x_curr(:, 2), 'c');
    plot(x_curr(:,3), x_curr(:, 4), 'r');
    plot(x_curr(:, 5), x_curr(:, 6), 'g');
end
xlabel('position')
ylabel('velocity')
axis square
% xl = [-Xmax, Xmax];
% yl = [-Xmax, Xmax];

% patch(R0*[-1, -1, 1, 1, -1], R0*[-1, 1, 1, -1, -1], 'm', 'EdgeColor', 'none');

% fcontour(@(x, y) 1 - cos(x).*cos(y) - sin(x).*sin(y), [xl, yl], 'levellist', [synclim], 'linecolor', 'r')
% xlim(xl);
% ylim(yl);



% title('Sampled Trajectories', 'FontSize', 16)
% 
figure(5)
clf
hold on
for i = 1:Ntraj
    x_curr = x_traj{i};
    plot(x_curr(:, 3) - x_curr(:, 1), x_curr(:, 5) - x_curr(:, 3), 'c');
%     % plot((cos(x_curr(:, 1))- cos(x_curr(:, 2))).^2 + (sin(x_curr(:, 1))- sin(x_curr(:, 2))).^2, 'c');
% 
% plot(1- (x_curr(:, 3).*x_curr(:, 4) + x_curr(:, 5).*x_curr(:, 6)), 'c');
% 
end

xlabel('$x_2 - x_1$', 'interpreter', 'latex')
ylabel('$x_3 - x_2$', 'interpreter', 'latex')

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
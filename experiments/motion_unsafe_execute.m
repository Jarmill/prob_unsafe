%2d motion when drift term is constant

%% variables and dynamics
t = sdpvar(1,1);
x = sdpvar(2,1);



% b = -0.1;
sigma = 0.1;
f =  [x(2); -(x(1) +x(2) + 0.5*x(1)^3)];
g = sigma * [0;1];
% f =  [-x(2); x(1)];
% f = [0; 0];
% f = 2*[0; -1];
% g = 2*[0; 0.1];

order_list = 4;
Norder = length(order_list);
v_list = cell(Norder, 3);
v0_list = zeros(Norder, 3);
obj_list = zeros(Norder, 3);

%% unsafe set

% theta_c = 5*pi/4; 
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

% x0 = [0; 0.75];
% x0 = [0.75; 0];
x0 = [0.85; -0.75]; %order 6 prb <= 3.0565354310e-01 (need to simulate)



%% Support Sets
% T = 1;
% T = 3;
% T = 5;
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
% order =6; 
for j = 1:length(order_list)
    order = order_list(j);
    d = 2*order;
    
    [v, cv, mv] = polynomial([t;x], d);
    [w, cw, mw] = polynomial([x], d);
    
    fT = T*f;
    gT = T*g;
    
    Lv = jacobian(v, t) + jacobian(v, x)*fT + 0.5*(gT)'*hessian(v, x)*(gT);
    
    v0 = replace(v, [t], [0]);
    
    %toggles for different experiments
    LEBESGUE = 1;
    CIRC = 0;
    
    consinit = [];
    coeffinit = [];
    
    v0x = replace(v0, x, x0);
    
    %leb
        leb = LebesgueBoxMom( d, [-Xmax, Xmax; -Xmax, Xmax]',1);
        w_leb = cw'*leb;
        objective_leb = w_leb;
    
    %circ
        gamma = sdpvar(1, 1);
            objective_circ = gamma;
            R0 = 0.2;
            X0 = struct('ineq', R0^2 - sum((x-x0).^2), 'eq', []);
            [put_init, consinit, coeffinit] = constraint_psatz(gamma-v0, X0, x, d);
    
    %point
            v0x = replace(v0, x, x0);
            objective_pt = v0x;
            
        % end
    % end
    
    zero_con = (coefficients(v0 - w, x)==0);
    
    [put_lie, conslie, coefflie] = constraint_psatz(-Lv, Xall, [t;x], d);
    [put_base, consbase, coeffbase] = constraint_psatz(v, Xall, [t;x], d);
    [put_unsafe, consunsafe, coeffunsafe] = constraint_psatz(v-1, Xuall, [t;x], d);
    
    cons = [conslie:'lie'; consbase:'base'; consunsafe:'unsafe'; zero_con:'zero'; consinit];
    coeff = [coefflie; coeffbase; coeffunsafe; cv; cw; coeffinit];
    
    opts = sdpsettings('solver', 'mosek');
    
    %% execute the routines
    %point
    [sol,u,Q] = solvesos(cons,objective_pt,opts,coeff);
    
    fprintf('point: %0.4d \n', value(objective_pt))
    % recovery
    v_rec_pt = value(cv)'*mv;
    v0_rec_pt = replace(v_rec_pt, t, 0);
    obj_rec_pt = value(objective_pt);
    
    
    %circ
    [sol,u,Q] = solvesos([cons; consinit],objective_circ,opts,[coeff; coeffinit]);
    
    % disp(value(objective))
    fprintf('circ: %0.4d \n', value(objective_circ))
    % recovery
    v_rec_circ= value(cv)'*mv;
    v0_rec_circ = replace(v_rec_circ, t, 0);
    obj_rec_circ = value(objective_circ);
    
    %circ
    [sol,u,Q] = solvesos([cons],objective_leb,opts,[coeff]);
    
    fprintf('leb: %0.4d \n', value(objective_leb))
    % recovery
    v_rec_leb= value(cv)'*mv;
    v0_rec_leb = replace(v_rec_leb, t, 0);
    obj_rec_leb = value(objective_leb);



    %% exports
    v_list{j, 1} = v_rec_pt;
    v_list{j, 2} = v_rec_circ;
    v_list{j, 3} = v_rec_leb;

    v0_list(j, :) = [v0_rec_pt, v0_rec_circ, v0_rec_leb];
    obj_list(j, :) = [obj_rec_pt, obj_rec_circ, obj_rec_leb];

    save("motion_unsafe_report.mat", "v0_list", "obj_list", "order_list");

end

using JuMP
using Mosek
using Random
using rational_peak_est
using Revise
using TSSOS
using SemialgebraicSets
using DynamicPolynomials


INIT_BOX = true;

#variables
@polyvar w[1:2];
@polyvar c[1:2];
@polyvar s[1:2];
@polyvar t

varsbase = [w; c; s];

#support sets

Wmax = 1;
X = Wmax.^2 .- w.^2;
circ = c.^2 .+ s.^2 .- 1;

# sigma = 0.15;
sigma = 0.05;
#parameter set 

#rational lift
# f = [-mu1*x[1] + alpha1*y[1]; -mu2*x[2] + alpha2*y[2]];
sin_diff = s[1]c[2] - c[1]*s[2];
fw1 = -s[1] -0.4*w[1] -0.5*(sin_diff);
fw2 = -0.5*s[2] - 0.5*(-sin_diff) - 0.5*w[2] + 0.05;
fw = [fw1; fw2];
f = [fw; s.*w; (-c).*w];
g = [1; 1; 0; 0; 0; 0] * sigma;


model = Model(optimizer_with_attributes(Mosek.Optimizer));

# set_optimizer_attribute(model, MOI.Silent(), true);

#problem
# order = 4;
order = 3;
gam =  @variable(model);
v, vc, vb = add_poly!(model, [varsbase; t], 2*order);

v = rem(v, circ);
# Tmax = 8;
Tmax = 5;




XT = [X; t*(1-t)];

#initial set
# R0 = 0.15;
# R0 = 0.05;
R0 = 0.01;
X0 = [c .- cos(R0); s.+ sin(R0); sin(R0) .- s];


#objective to maximize
p = 1 - c[1]*c[2] - s[1]*s[2];

synclim = 0.2;
# synclim = 0.9;
XuT = [X; p - synclim; t*(1-t)];
# XuT = [X;  t*(1-t)];

#Lie derivative
# deg_change = 1; #floor(degree f/ 2), but I don't know how to take the degree of the polynomial array f
deg_change = 1;
# model,info1 = add_psatz!(model, -Lv, [x; t], XT, [], d, QUIET=true, CS=true, TS="block", Groebnerbasis=false)
dvdx = differentiate(v, varsbase);
fterm = sum(f .* dvdx);
v_hess = differentiate.(dvdx', varsbase)
gterm = g'*v_hess*g;
Lv = differentiate(v, t) + Tmax * fterm + (Tmax/2)*gterm;
Lv = rem(Lv, circ);
# model,info1 = add_psatz!(model, -Lv, [varsbase; t], [XT; synclim-p], circ, order+deg_change, QUIET=true, CS=true, TS="block", Groebnerbasis=true);
model,info1 = add_psatz!(model, -Lv, [varsbase; t], XT, circ, order+deg_change, QUIET=true, CS=true, TS="block", Groebnerbasis=true);

model,info2 = add_psatz!(model, v, [varsbase; t], [XT; synclim - p], circ, order, QUIET=true, CS=true, TS="block", Groebnerbasis=true);
# model,info2 = add_psatz!(model, v - p, [varsbase; t], [XT; synclim - p], circ, order, QUIET=true, CS=true, TS="block", Groebnerbasis=true);
# model,info2 = add_psatz!(model, v - p, [varsbase; t], XT, circ, order, QUIET=true, CS=true, TS="block", Groebnerbasis=true);


GB = groebnerbasis([circ;p-synclim]);
v_u = rem(v, GB);
# model,info3 = add_psatz!(model, v-1, [varsbase; t], [XT], GB, order, QUIET=true, CS=true, TS="block", Groebnerbasis=true);
# model,info3 = add_psatz!(model, v-1, [varsbase; t], XuT, circ, order, QUIET=true, CS=true, TS="block", Groebnerbasis=true);
model,info3 = add_psatz!(model, v_u-1, [varsbase; t], XT, GB, order, QUIET=true, CS=true, TS="block", Groebnerbasis=true);

# #initial condition
if ~INIT_BOX
    v0 = subs(v, [t; w; c; s]=>[0; 0; 0; 1; 1; 0; 0]); #subs(p, y=>x^2)
    @constraint(model, gam-v0>=0);
else    
    v0 = subs(v, [t; w]=>[0; 0; 0]); #subs(p, y=>x^2)
    model,info4 = add_psatz!(model, gam - v0, [c; s], [X0; synclim - p], circ, order, QUIET=true, CS=true, TS="block", Groebnerbasis=true);
end
@objective(model, Min, gam)
optimize!(model)
status = termination_status(model)
if status != MOI.OPTIMAL
    println("termination status: $status")
    status = primal_status(model)
    println("solution status: $status")
end
objv = objective_value(model)
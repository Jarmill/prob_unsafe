
#Author: Jie Wang, https://github.com/wangjie212/SparseDynamicSystem/blob/main/example/example.jl
using JuMP
using Revise
using MosekTools
using DynamicPolynomials
using MultivariatePolynomials
using SparseDynamicSystem
#variables
n = 2
@polyvar t
@polyvar x[1:n]

#system
# f = [x[2]; -x[1] - x[2] + x[1]^3/3]
# f = [x[2]; (-x[1] - x[2] - x[1]^3)/2];
f = [x[2]; -x[1] - x[2] - 0.5*x[1]^3];
sigma = 0.1;
g = [0; sigma];
# f = -x;

#support sets
Tmax = 5;

x0 = [0.85; -0.75];

Box = 2.5;
X = Box^2 .- [x[1]^2; x[2]^2];
XT = [X; t*(1-t)]
Xu = [-0.5625-x[1]-1.5*x[2]-x[1]^2-x[2]^2; -0.75-x[2]];
XuT = [Xu; X; t*(1-t)];
# X0 = [R0^2 - sum(x-C0).^2)]


#order of polynomials/moments (2d)
d = 6; 

# find the probability of unsafety


#model and auxiliary functions
model = Model(optimizer_with_attributes(Mosek.Optimizer))
set_optimizer_attribute(model, MOI.Silent(), true)
v, vc, vb = add_poly!(model, [t; x], 2d)
gamma = @variable(model); 

#Lie derivative
deg_change = 1; #floor(degree f/ 2), but I don't know how to take the degree of the polynomial array f
# deg_change = 0;
vx = differentiate(v, x);
vxx = differentiate(vx, x);

fT = f*Tmax;
gT = g*Tmax;

Lv = differentiate(v, t) + sum(fT .* vx) + ((vxx*gT)')*gT/2; 
model,info1 = add_psatz!(model, -Lv, [t; x], XT, [], d+deg_change, QUIET=true, CS=true, TS="block", Groebnerbasis=false)

#v above the indicator function of Xu
model,info2 = add_psatz!(model, v, [t; x], XT, [], d, QUIET=true, CS=true, TS="block", Groebnerbasis=false)
model,info3 = add_psatz!(model, v-1, [t; x], XuT, [], d, QUIET=true, CS=true, TS="block", Groebnerbasis=false)

#initial condition
# v0 = subs(v, [t; x]=>[0; x0]); #subs(p, y=>x^2)
v0 = v(t=>0, x=>x0); #subs(p, y=>x^2)

@objective(model, Min, v0)
optimize!(model)
status = termination_status(model)
if status != MOI.OPTIMAL
    println("termination status: $status")
    status = primal_status(model)
    println("solution status: $status")
end
objv = objective_value(model)

v_rec = value.(vc)'*vb;
print(objv)
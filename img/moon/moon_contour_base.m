function [x_cont] = moon_contour_base(h_in, h_out, dist, Narc)
%MOON_CONTOUR_BASE the L2 contour of a moon.
%The moon contains points [-1, 0] and [1, 0]. It is defined between a
%circle at height h_in and at height h_out, where the height is the lowest
%y coordinate on the circle.
%
% x_moon is the set of points with an L2 distance of 'dist' away from the
% moon. This is a contour.

%Input:
%   h_in:   height of inner circle
%   h_out:  height of outer circle
%   dist:   distance of contour
%   Narc:   number of sample pointsof arc
%
%Output:
%   x_cont: contour of the moon
%% input processing
assert(h_out > h_in, 'Outer circle height should be larger than inner circle height');

if nargin < 4
    Narc = 201;
end


%conditions:
% nominal:  h_in > 0, dist <= h+c
% pinch:    h_in > 0, dist >  h+c
% flat:     h_in = 0

%centers and radii of the circle
c_in = 0.5*(1/h_in - h_in);
r_in = 0.5*(1/h_in + h_in);

c_out = 0.5*(1/h_out - h_out);
r_out = 0.5*(1/h_out + h_out);

%% intersection computation
syms x y

eq_in = ((r_in-dist)^2 == x^2 + (y-c_in)^2);
eq_side = ((dist)^2 == (x-1)^2 + (y)^2);
eq_out = ((r_out+dist)^2 == x^2 + (y-c_out)^2);

%intersection between side and outer
[xout, yout] = solve([eq_out,eq_side], [x,y]);
if h_in == 0
    %degenerate: flat top of moon
    xin = 1;
    yin = dist;
    theta_inner = NaN;
else
    %curved inner circle
    if dist > h_in+c_in
        %pinch point at x=0
        xin = 0;
        yin = sqrt(dist^2 - 1);
        theta_inner = NaN;
    else
        %nominal
        [xin, yin] = solve([eq_in,eq_side], [x,y]);
        
        xin = real(xin(1));
        yin = real(yin(1));
        theta_inner = atan2(yin-c_in, xin);
    end           
end
    
%% figure out angle ranges
theta_outer = atan2(yout-c_out, xout);
theta_outer_side = atan2(yout, xout-1);
theta_inner_side = atan2(yin, xin-1);

Nout = floor(Narc * (theta_outer+pi/2)/(2*pi))+1;
Nside = floor(Narc * (theta_outer_side - theta_inner_side)/(2*pi))+1;

theta_out = linspace(-pi/2, theta_outer, Nout);
theta_side = linspace(theta_outer_side,wrapTo2Pi(theta_inner_side), Nout);

x_out = (r_out+dist)*[cos(theta_out); sin(theta_out)] + [0; c_out];
x_side = (dist)*[cos(theta_side); sin(theta_side)] + [1; 0];


if isnan(theta_inner)
    x_in = [];
else
    Nin = floor(Narc * (theta_inner+pi/2)/(2*pi))+1;
    theta_in = linspace(theta_inner, -pi/2, Nout);
    x_in = (r_in-dist)*[cos(theta_in); sin(theta_in)] + [0; c_in];
end


%% assemble the moon contour

x_half = [x_out, x_side, x_in];

x_refl = diag([-1,1])*x_half(:, end:-1:1);

x_cont = double([x_half, x_refl]);

end


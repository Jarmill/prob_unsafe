function [x_cont] = moon_contour_base_2(h_in, h_out, dist, Narc)
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
%   Narc:   number of sample points of arc
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

%intersection between side and outer
%   x_inner:    intersection between inner and side circle
%   x_outer:    intersection between outer and side circle 
[x_inner, x_outer] = moon_intersection_points(c_in, c_out, dist);



if h_in == 0
    theta_inner = NaN;
else
    %curved inner circle
    if dist > h_in+c_in
        %pinch point at x=0
        x_inner(1) = 0;
        x_inner(2) = sqrt(dist^2 - 1);
        theta_inner = NaN;
    else
        %nominal
        theta_inner = atan2(x_inner(2)-c_in, x_inner(1));
    end           
end
    
%% figure out angle ranges
theta_outer = atan2(x_outer(2)-c_out, x_outer(1));
%outer circle arc
Nout = floor(Narc * (theta_outer+pi/2)/(2*pi))+1;
theta_out = linspace(-pi/2, theta_outer, Nout);
x_out = (r_out+dist)*[cos(theta_out); sin(theta_out)] + [0; c_out];

%inner circle arc
if isnan(theta_inner)
    x_in = [];
else
    Nin = floor(Narc * (theta_inner+pi/2)/(2*pi))+1;
    theta_in = linspace(theta_inner, -pi/2, Nout);
    x_in = (r_in-dist)*[cos(theta_in); sin(theta_in)] + [0; c_in];
end

%side circle arc
theta_outer_side = atan2(x_outer(2), x_outer(1)-1);
theta_inner_side = atan2(x_inner(2), x_inner(1)-1);


if (dist > 1) && (h_in > 1)
    sdist = sqrt(dist^2-1);
    theta_side_side_top = wrapTo2Pi(atan2(sdist, -1));
    theta_side_side_bot = wrapTo2Pi(atan2(-sdist, -1));
    
    
    
    Nside_top = floor(Narc * (theta_side_side_top-theta_outer_side)/(2*pi))+1;
    theta_side_top = linspace(theta_outer_side,(theta_side_side_top), Nside_top);
    x_side_top = (dist)*[cos(theta_side_top); sin(theta_side_top)] + [1; 0];
    %the side circles overlap. There may be 2 connected components
    
    x_side = x_side_top;
    if ~isnan(theta_inner) % (dist < r_in)
            Nside_bot = floor(Narc * (theta_side_side_bot - theta_inner_side)/(2*pi))+1;
            theta_side_bot = linspace(theta_side_side_bot,wrapTo2Pi(theta_inner_side), Nside_bot);
            x_side_bot = (dist)*[cos(theta_side_bot); sin(theta_side_bot)] + [1; 0];
            
            
        x_side = [x_side, [NaN; NaN], x_side_bot];
    end    
else
    %the side circles do not overlap, at most 1 connected component
    Nside = floor(Narc * (theta_outer_side - theta_inner_side)/(2*pi))+1;
    theta_side = linspace(theta_outer_side,wrapTo2Pi(theta_inner_side), Nout);
    x_side = (dist)*[cos(theta_side); sin(theta_side)] + [1; 0];
end


%% assemble the moon contour

x_half = [x_out, x_side, x_in];
% x_half = [x_out, [NaN; NaN], x_in];

x_refl = diag([-1,1])*x_half(:, end:-1:1);

x_cont = double([x_half, x_refl]);
% x_cont = x_half;

end


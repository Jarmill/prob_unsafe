function [x_inner, x_outer] = moon_intersection_points(c_in, c_out, dist)
%MOON_INTERSECTION_POINTS find the points of intersection between the
%inner/outer circles and the corner circles for the moon contour

%Input:
%   c_in:   center of inner circle
%   c_out:  center of outer circle
%   dist:   L2 distance to circle

%Output:
%   x_inner:    intersection between inner and side circle
%   x_outer:    intersection between outer and side circle

poly_center = @(c) [c^2+1, -2*(c^2+1), c^2 + 1 - dist^2];

%inner/side
if c_in == Inf
    x_in = 1;
    y_in = dist;
else
    poly_in = poly_center(c_in);
    roots_in = roots(poly_in);
    x_in = min(roots_in);
    y_in = c_in*(1-x_in);
end
%outer/side point
poly_out = poly_center(c_out);
roots_out = roots(poly_out);
x_out = max(roots_out);
y_out = c_out*(1-x_out);

x_inner = [x_in; y_in];
x_outer = [x_out; y_out];



end


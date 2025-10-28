function pts = get_ring_points(n, rad, center)
% n is an integer specifying number of points
% rad is a float
% center (optional) specifies the center of the ring
if nargin < 3
    center = [0 0];
end

ring_thetas = linspace(0, 2*pi, n);
ring_xvals = rad * cos(ring_thetas) + center(1);
ring_yvals = rad * sin(ring_thetas) + center(2);
pts = [ring_xvals; ring_yvals];
end

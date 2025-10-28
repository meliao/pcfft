function pts = get_sphere_points(n, radius, center)
    % Generates n points approximately uniformly distributed over the 
    % sphere. 
    % From: https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere

    if nargin < 2
        radius = 1.0;
    end
    if nargin < 3
        center = zeros(3, 1);
    end
    % Golden angle in radians
    phi = pi * (sqrt(5) - 1);

    idxes = 0:n-1;


    yvals = linspace(1 / n- 1, 1 - 1 / n, n);
    r = sqrt(1 - yvals .* yvals);
    disp(min(r));
    thetas = phi * idxes;
    xvals = cos(thetas) .* r;
    zvals = sin(thetas) .* r;

    pts = [xvals; yvals; zvals];
    pts = pts * radius + center;

end

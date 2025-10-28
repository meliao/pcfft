% addpath(genpath('../../pcfft'));
% addpath(genpath('~/FLAM'))

nsrc = 1e4;
L = 15;
xs = L*(2*rand(2,nsrc) - 1);
% xs = xs .* vecnorm(xs).^2;

figure(1);clf
scatter(xs(1,:), xs(2,:),'.')


crad = 2;

n_nbr = 1e2;
half_side = spread_halfside(xs, n_nbr, crad);
half_side

pts = half_side * [[-1;-1],[-1;1], [1;1], [1;-1] , [-1;-1]];

thetas = linspace(0,2*pi);
pxy2 = crad*[cos(thetas);sin(thetas)];

hold on
plot(pts(1,:), pts(2,:),'linewidth',2)
plot(pxy2(1,:), pxy2(2,:),'linewidth',2)
hold off

axis equal


rx = xs(1,:) - xs(1,:).';
ry = xs(2,:) - xs(2,:).';
r = sqrt(rx.^2 + ry.^2);
r = sort(r(:));

% check that n_nbr is about right
[sum(r < crad*2*half_side) / nsrc, n_nbr]

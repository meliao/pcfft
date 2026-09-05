addpath(genpath('../../pcfft'));
addpath(genpath('~/software/chunkie'));
addpath(genpath('~/software/fmm3dbie'));

crad = 2;
n_nbr_list = round(logspace(1, 3, 10));

cases = {};

cparams = [];
cparams.maxchunklen = 0.02;
pref = [];
pref.k = 16;
chnkr = chunkerfunc(@(t) starfish(t,5,0.3,[0;0],0,1), cparams, pref);
cases{end+1} = struct('dim',2,'name','oversampled starfish','xs',chnkr.r(:,:));

cparams2 = [];
cparams2.maxchunklen = 0.3;
pref2 = [];
pref2.k = 16;
L = 12;
ntry = 1000;
nscat = 30;
ctrs = zeros(2,1);
chnkrs = [];
for i = 1:nscat
    narms = randi([3,6]);
    amp = 0.2+0.2*rand();
    phi = 2*pi*rand();
    chnkr_i = chunkerfunc(@(t) starfish(t,narms,amp,ctrs(:,i),phi,0.5), cparams2, pref2);
    chnkrs = [chnkrs, chnkr_i];
    for j = 1:ntry
        theta = 2*pi*rand();
        tmp = L*rand()*[cos(theta);sin(theta)];
        rmin = min(vecnorm(tmp - ctrs));
        if rmin > 2.5; break; end
    end
    if j == ntry; error('Could not place next boundary'); end
    ctrs = [ctrs, tmp];
end
chnkrs = merge(chnkrs);
cases{end+1} = struct('dim',2,'name','many starfish','xs',chnkrs.r(:,:));

S = geometries.ellipsoid([1.5,1,0.7], [12,10,8], [0;0;0], 8, 11);
cases{end+1} = struct('dim',3,'name','oversampled ellipsoid','xs',S.r(:,:));

L3 = 10;
ntry = 1000;
nscat3 = 15;
ctrs3 = zeros(3,1);
surfs = [];
for i = 1:nscat3
    abc = 0.5+0.5*rand(1,3);
    surfi = geometries.ellipsoid(abc, [3,3,3], ctrs3(:,i), 6, 11);
    surfs = [surfs, surfi];
    for j = 1:ntry
        tmp = L3*(2*rand(3,1)-1);
        rmin = min(vecnorm(tmp - ctrs3));
        if rmin > 3; break; end
    end
    if j == ntry; error('Could not place next boundary'); end
    ctrs3 = [ctrs3, tmp];
end
surfs = merge(surfs);
cases{end+1} = struct('dim',3,'name','many ellipsoids','xs',surfs.r(:,:));

for c = 1:numel(cases)
    dim = cases{c}.dim;
    name = cases{c}.name;
    xs = cases{c}.xs;
    nsrc = size(xs,2);

    tree = createns(xs.');

    actual_nbrs = zeros(size(n_nbr_list));
    for k = 1:numel(n_nbr_list)
        n_nbr = n_nbr_list(k);
        half_side = spread_halfside(xs, n_nbr, crad);
        radius = 2*crad*half_side;

        idx = randperm(nsrc, min(2000,nsrc));
        counts = rangesearch(tree, xs(:,idx).', radius);
        nc = cellfun(@numel, counts) - 1;
        actual_nbrs(k) = mean(nc);
        fprintf('%s (dim=%d, npts=%d): n_nbr=%d half_side=%.4f radius=%.4f actual=%.2f\n', ...
            name, dim, nsrc, n_nbr, half_side, radius, actual_nbrs(k));
    end

    figure(c); clf
    loglog(n_nbr_list, actual_nbrs, 'o-', 'linewidth', 2)
    hold on
    loglog(n_nbr_list, n_nbr_list, 'k--')
    hold off
    xlabel('n_{nbr}')
    ylabel('actual mean neighbors')
    title(sprintf('%s (dim=%d, npts=%d)', name, dim, nsrc))
    legend('measured', 'y = x', 'location', 'northwest')
    grid on
end

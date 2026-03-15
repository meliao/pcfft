% Demonstrate application of the PCFFT method to a Gaussian kernel
%% setup problem
nsrcs = 2e4;
ntargs = 2e4;
dim = 2;

L = 2;

% Define the sources
srcs = [];
srcs.r = L*2*(rand(dim,nsrcs) - 0.5);


% Define the targets
targs = [];
targs.r = L*2*(rand(dim,ntargs) - 0.5);

covar = [2,1;1,3];

kern_0 = @(s,t) Gaussian_kern(s, t, covar);

% charges
str = randn(nsrcs,1);

%% Precomputation

eps = 1e-7;

t1 = tic;

% this is not a PDE kernel so we need proxy against values at a few radii
proxy_opts = [];
proxy_opts.multi_shells = true;

% setup grid
[grid_info, proxy_info] = get_grid(kern_0, srcs, targs, eps, [],proxy_opts);
% get spreading operators
[A_spread_s, sort_info_s]= get_spread(kern_0, [], srcs, ...
    grid_info, proxy_info);
[A_spread_t, sort_info_t]= get_spread(kern_0, [], targs, ...
    grid_info, proxy_info);
% build corrections
A_addsub = get_addsub(kern_0, [], srcs, targs, ...
    grid_info, proxy_info, sort_info_s, sort_info_t, A_spread_s, A_spread_t);
% get DFT of kernel
kern_0hat = get_kernhat(kern_0,grid_info);

tprecom = toc(t1);

%% True solution
tic;
Atrue = kern_0(srcs,targs);
t_densecomp = toc;
tic;
utrue = Atrue*str;
tdens_app = toc;

%% Apply


tic;
u = pcfft_apply(str,A_spread_s,A_spread_t,A_addsub,kern_0hat);
tapply = toc;
err = norm(u - utrue) / norm(utrue);
fprintf('compression error = %e, tol = %.2e\n', err,eps)
fprintf('pcfft precomputation time = %.2e, apply time = %.2e\n', tprecom,tapply)
fprintf('Dense precomputation time = %.2e, apply time = %.2e\n', t_densecomp,tdens_app)

%% Plot kernel
xx = linspace(-L,L);
[X,Y] = meshgrid(xx,xx);
targs = []; targs.r = [X(:).';Y(:).'];
src = []; src.r = [0;0];
u = kern_0(src,targs);
figure(1);
h = pcolor(X,Y,reshape(u,size(X))); h.EdgeColor = 'none'; colorbar
title('$K(0,\cdot)$','Interpreter','latex')
set(gca,'FontSize',16)

function val = Gaussian_kern(src_pts,target_pts,covar)
% src_pts has shape (2, M)
% target_pts has shape (2, N)
% Computes log{|| src - target||}
% Output shape is (N, M)

% Shape (N, M)
dist = 0;
for i = 1:size(src_pts.r(:,:),1)
    for j = 1:size(src_pts.r(:,:),1)
    dist = dist + covar(i,j)*(src_pts.r(i, :) - target_pts.r(i, :).').*(src_pts.r(j, :) - target_pts.r(j, :).');
    end
end

val = exp(-dist/2);
end
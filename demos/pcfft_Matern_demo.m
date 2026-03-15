% Demonstrate application of the PCFFT method to the Matern kernel

addpath(genpath("../pcfft"));

%% setup problem
nsrcs = 2e4;
ntargs = 2e4;
dim = 3;

L = 2;

% Define the sources
srcs = [];
srcs.r = L*2*(rand(dim,nsrcs) - 0.5);


% Define the targets
targs = [];
targs.r = L*2*(rand(dim,ntargs) - 0.5);

nu = 1.5;
kern_0 = @(s,t) Matern_kern(s, t, nu);

% charges
str = randn(nsrcs,1);

%% Precomputation

eps = 1e-6;

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

function val = Matern_kern(src_pts,target_pts,nu)
% src_pts has shape (2, M)
% target_pts has shape (2, N)
% Computes log{|| src - target||}
% Output shape is (N, M)

% Shape (N, M)
dist = 0;
for i = 1:size(src_pts.r(:,:),1)
    dist = dist + (src_pts.r(i, :) - target_pts.r(i, :).').^2;
end

dist = sqrt(dist);

val = (dist).^nu.*besselk(nu,dist);
end
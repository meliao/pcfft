%% setup problem
nsrcs = 2e4;
ntargs = 2e4;

L = 2;

srcs = [];
srcs.r = L*2*(rand(2,nsrcs) - 0.5);
srcs.n = randn(2,nsrcs);
srcs.n = srcs.r;
targs = [];
targs.r = L*2*(rand(2,ntargs) - 0.5);
targs.n = randn(2,ntargs);

kern_0 = @(s,t) log_kernel(s, t);
% kern_0 = @(s,t) one_over_r_kernel2D(s, t);
kern_s = @(s,t) wrap_d(kern_0,s,t);
% kern_s = kern_0;
kern_t = kern_0;
kern_st = kern_s;

% charges
str = randn(nsrcs,1);

%% Precompuation

eps = 1e-6;

t1 = tic;

% setup grid
[grid_info, proxy_info] = get_grid(kern_0, srcs, targs, eps);
% get spreading operators
[A_spread_s, ~, sort_info_s]= get_spread(kern_0, kern_s, srcs, ...
    grid_info, proxy_info, {'r','n'});
[A_spread_t, ~, sort_info_t]= get_spread(kern_0, kern_t, targs, ...
    grid_info, proxy_info);
% build corrections
[A_addsub] = get_addsub(kern_0, kern_st, srcs, targs, ...
    grid_info, proxy_info, sort_info_s, sort_info_t, A_spread_s, A_spread_t);
% get DFT of kernel
kern_0hat = get_kernhat(kern_0,grid_info);

tprecom = toc(t1)

%% True solution
tic;
Atrue = kern_st(srcs,targs);
true = toc
tic;
utrue = Atrue*str;
tdens_app = toc

%% Apply


tic;
u = pcfft_apply(str,A_spread_s,A_spread_t,A_addsub,kern_0hat);
tapply = toc
err = norm(u - utrue) / norm(utrue)

function val = wrap_d(kern0,s,t)

[~,grad] = kern0(s,t);

val = grad(:,:,1).*s.n(1,:) + grad(:,:,2).*s.n(2,:);

end
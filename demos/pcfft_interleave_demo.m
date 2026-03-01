% In this demo we evaluate 
%     u(x) = sum_j K(x,y_j) mu_j + d_j . grad_y K(x,y_j) rho_j
% and its gradient for K(x,y) = log(||x-y||)

%% setup problem
nsrcs = 2e3;
ntargs = 2e3;

L = 2;

% Define the sources
srcs = [];
srcs.r = L*2*(rand(2,nsrcs) - 0.5);
% Specify the normal vector at each source point. 
srcs.n = randn(2,nsrcs);
srcs.d = randn(2,nsrcs);
srcs = chnkr;
nsrcs = chnkr.npt;

% Define the targets
targs = [];
targs.r = L*2*(rand(2,ntargs) - 0.5);
targs.n = randn(2,ntargs);
targs.d = randn(2,ntargs);
targs = chnkr;

kern_0 = @(s,t) log_kernel(s, t);
% Our implementation of log_kernel can return the gradient; see def of 
% wrap_d at the bottom of this file.
kern_s = kern_0;
kern_t = @(s,t) my_kernel(t,s).';
kern_st = @(s,t) my_kernel(s,t);

zk = 0;
skern = @(s,t) chnk.flex2d.kern(zk, s, t, 's');   
ekern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');   
bkern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_bcs');   
fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate');

kern_0 = skern;
kern_s = skern;
kern_t = @(s,t) bkern(t,s).';
kern_st = bkern;

kern_0 = skern;
kern_s = ekern;
kern_t = skern;
kern_st = ekern;

kern_0 = skern;
kern_s = ekern;
kern_t = @(s,t) bkern(t,s).';
kern_st = fkern;

% ekern2 = @(src,targ) kern_component(ekern, src, targ,1,1,[1,2]);
% kern_0 = skern;
% kern_s = ekern2;
% kern_t = skern;
% kern_st = ekern2;

% bkern2 = @(src,targ) kern_component(bkern, src, targ,2,1,[2,1]);
% kern_0 = skern;
% kern_s = skern;
% kern_t = @(s,t) bkern2(t,s).';
% kern_st = bkern2;

% charges
str = randn(nsrcs,1);
str = randn(2*nsrcs,1);

%% Precompuation

eps = 1e-7;

t1 = tic;

% setup grid
[grid_info, proxy_info] = get_grid(kern_0, srcs, targs, eps);
% get spreading operators
[A_spread_s, sort_info_s]= get_spread(kern_0, kern_s, srcs, ...
    grid_info, proxy_info, {'r','n','d'});
[A_spread_t, sort_info_t]= get_spread(kern_0, kern_t, targs, ...
    grid_info, proxy_info, {'r','n','d'});
% build corrections
A_addsub = get_addsub(kern_0, kern_st, srcs, targs, ...
    grid_info, proxy_info, sort_info_s, sort_info_t, A_spread_s, A_spread_t);
% get DFT of kernel
kern_0hat = get_kernhat(kern_0,grid_info);

tprecom = toc(t1)
%% True solution
tic;
Atrue = kern_st(srcs,targs);
Atrue(1:2*size(Atrue,1)+2:end) = 0;
Atrue(2:2*size(Atrue,1)+2:end) = 0;
Atrue(1+size(Atrue,1):2*size(Atrue,1)+2:end) = 0;
Atrue(2+size(Atrue,1):2*size(Atrue,1)+2:end) = 0;
true = toc
tic;
utrue = Atrue*str;
tdens_app = toc

%% Apply


tic;
u = pcfft_apply(str,A_spread_s,A_spread_t,A_addsub,kern_0hat);
tapply = toc
err = norm(u - utrue) / norm(utrue)

% function val = log_kern(src,targ)
% 
% 
% 
% if strcmp(type, 's')
% 
% 
% elseif strcmp(type, 'd')
% 
% elseif strcmp(type, 'sgrad')
% 
% elseif strcmp(type, 'dgrad')
% 
% 
% end
% end



function kvals = log_kernel(src_pts, target_pts)
% src_pts has shape (2, M)
% target_pts has shape (2, N)
% Computes log{|| src - target||}
% Output shape is (N, M)

% Shape (N, M)
rx = src_pts.r(1, :) - target_pts.r(1, :).';
ry = src_pts.r(2, :) - target_pts.r(2, :).';

dist = sqrt(rx.^2 + ry.^2);

kvals = log(dist);

end

function vals = my_kernel(src_pts, target_pts)
% src_pts has shape (2, M)
% target_pts has shape (2, N)
% Computes grad log{|| src - target||}
% Output shape is (N, M)

% Shape (N, M)
rx = src_pts.r(1, :) - target_pts.r(1, :).';
ry = src_pts.r(2, :) - target_pts.r(2, :).';

dist = sqrt(rx.^2 + ry.^2);

kvals = log(dist);
grad = zeros(2,size(rx,1),size(rx,2));
grad(1,:,:) = -rx ./ dist.^2;
grad(2,:,:) = -ry ./ dist.^2;


% hess = zeros(3,size(rx,1),size(rx,2));
% 
% hess(1,:,:) = 2*rx2./(r4)-1.0./(r2);
% hess(2,:,:) = 2*rx.*ry./(r4);
% hess(3,:,:) = 2*ry2./(r4)-1.0./(r2);
% 
% 
% vals1 = [reshape(kvals,1,size(rx,1),size(rx,2)); reshape(grad,2,size(rx,1),size(rx,2))];
% vals1 = reshape(vals1, [], size(rx,2));


vals = [reshape(kvals,1,size(rx,1),size(rx,2)); reshape(grad,2,size(rx,1),size(rx,2))];
vals = reshape(vals, [], size(rx,2));

% vals = kvals;
end
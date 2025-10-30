addpath(genpath('../../pcfft'));

rad = 10.0;
tol = 1e-13;
dim = 2;
boxhalf_sidelen = 0.5;

src_info_2d = struct;
n_src = 13;
src_info_2d.r = [2;1].*(rand(2, n_src) - 0.5) * boxhalf_sidelen;
src_info_2d.weights = rand(n_src, 1);

targ_info_2d = struct;
targ_info_2d.radius = 4.0;
n_targ = 17;
targ_info_2d.r = (rand(2, n_targ) - 0.5) * boxhalf_sidelen;

[grid_info, proxy_info] = get_grid(@log_kernel,src_info_2d, targ_info_2d, tol);

rgrid = grid_info.r;
ngrid = grid_info.ngrid;
Lbd = grid_info.Lbd;

% define strengths and direct u

str = zeros(1,size(rgrid,2));

% ipt = ngrid(2)+1 + (2*ngrid(2)+1)*(ngrid(1)+1);
ipt = 1;
str(ipt) = 1;
kerns = log_kernel(rgrid(:,ipt), rgrid);

u = kerns*str(ipt);

% compute u using fft
kern_hat = get_kernhat(@log_kernel,rgrid,ngrid, Lbd, grid_info.dx);

str_hat = fftn(reshape(str,size(kern_hat)));
u_hat = kern_hat .* str_hat;
ugrid = ifftn(u_hat);

%compare
ugrid = ugrid(:);

% get all 'real' grid points
[idx, idy] = meshgrid(1:ngrid(1), 1:ngrid(2));
idreal = idy(:) + (2*ngrid(2)+1)*idx(:);

assert(norm(u(idreal) - ugrid(idreal))<1e-10)

%%

figure(1)
subplot(1,3,1)
scatter(rgrid(1,:), rgrid(2,:), [], u(:))
colorbar
subplot(1,3,2)
scatter(rgrid(1,:), rgrid(2,:), [], ugrid(:))
colorbar
subplot(1,3,3)
scatter(rgrid(1,:), rgrid(2,:), [], log10(abs(ugrid(:)-u(:))))
colorbar

% rpad = 2;
% xx = Lbd(1,1) + (0:rpad*ngrid(1)) *grid_info.dx;
% yy = Lbd(2,1) + (0:rpad*ngrid(2)) *grid_info.dx;
% [X, Y] = meshgrid(xx,yy);
% 
% figure(4);
% h = pcolor(X,Y, reshape(u,size(X))); h.EdgeColor = 'None';



%% 3d
dim = 3;
src_info_3d = struct;
n_src = 13;
src_info_3d.r = (rand(3, n_src) - 0.5) * half_sidelen;
src_info_3d.weights = rand(n_src, 1);

targ_info_3d = struct;
targ_info_3d.radius = 4.0;
targ_info_3d.r = (rand(3, n_targ) - 0.5) * half_sidelen;
tol = 1e-08;


[grid_info, proxy_info] = get_grid(@one_over_r_kernel,src_info_3d, targ_info_3d, tol);


rgrid = grid_info.r;
ngrid = grid_info.ngrid;
Lbd = grid_info.Lbd;
% define strengths and direct u

str = zeros(1,size(rgrid,2));

% ipt = ngrid(2)+1 + (2*ngrid(2)+1)*(ngrid(1)+1);
ipt = 1;
str(ipt) = 1;
kerns = one_over_r_kernel(rgrid(:,ipt), rgrid);

u = kerns*str(ipt);
% %%
% compute u using fft
kern_hat = get_kernhat(@one_over_r_kernel,rgrid,ngrid, Lbd, grid_info.dx);

str_hat = fftn(reshape(str,size(kern_hat)));
u_hat = kern_hat .* str_hat;
ugrid = ifftn(u_hat);

%compare
ugrid = ugrid(:);

% get all 'real' grid points
[idx, idy, idz] = meshgrid(1:ngrid(1), 1:ngrid(2), 1:ngrid(3));
idreal = idz(:) + (2*ngrid(3)+1)*(idy(:) + (2*ngrid(2)+1)*idx(:));
% idreal = 1:21;
norm(u(idreal) - ugrid(idreal))

% %%
% 
% ngrid = 4e2 * [1;1;1]';
% % ngrid = 512 * [1;1;1]';
% 
% u = randn(ngrid);
% tic;
% uhat = fftn(u, 2*ngrid);
% toc;
% 
% u2 = randn(2*ngrid);
% tic;
% uhat2 = fftn(u2);
% toc;

%%


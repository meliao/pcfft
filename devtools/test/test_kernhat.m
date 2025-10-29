addpath(genpath('../../pcfft'));

rad = 10.0;
tol = 1e-13;
dim = 2;
half_sidelen = 0.5;

k = @(s,t) log_kernel(s,t);

src_info_2d = struct;
n_src = 13;
src_info_2d.r = [2;1].*(rand(2, n_src) - 0.5) * half_sidelen;
src_info_2d.weights = rand(n_src, 1);

targ_info_2d = struct;
targ_info_2d.radius = 4.0;
n_targ = 17;
targ_info_2d.r = (rand(2, n_targ) - 0.5) * half_sidelen;


[grid_info, proxy_info, rgrid, ngrid, Lbd] = get_grid(@log_kernel,src_info_2d, targ_info_2d, tol);



str = zeros(1,size(rgrid,2));

ipt = 80 + 181*80;
str(ipt) = 1;
kerns = log_kernel(rgrid(:,ipt), rgrid);

u = kerns*str(ipt);

kern_hat = get_kernhat(@log_kernel,rgrid,ngrid, Lbd, grid_info.dx);


str_hat = fft2(reshape(str,size(kern_hat)));
u_hat = kern_hat .* str_hat;
ugrid = ifft2(u_hat);

norm(u(:) - ugrid(:))

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



figure(2)
scatter(rgrid(1,:), rgrid(2,:),'.')




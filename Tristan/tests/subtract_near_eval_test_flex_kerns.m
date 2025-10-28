addpath(genpath('~/software/chunkie'))
addpath(genpath('~/projects/fmm3dbie/matlab'))

% norder = 8;
% rad = 15;
% S = geometries.disk([rad,rad],[],[3 3 6],norder);
% nch = 4*6;
% cparams.ta = pi/nch;
% cparams.tb = 2*pi + cparams.ta;
% chnkr = chunkerfuncuni(@(t) ellipse(t),nch,cparams);
% chnkr = chnkr.move([0;0],[0;0],0,rad);
% kappa = signed_curvature(chnkr);
% kappa = kappa(:);

rad = 2;
norder = 8;
S = geometries.disk([rad,rad],[],[6 4 8],norder); % [6 4 8] 
chnkr = chunkerfuncuni(@(t) ellipse(t),32);
chnkr = chnkr.move([0;0],[0;0],pi/32,2);
chnkr = sort(chnkr);

chnkr.rstor = [chnkr.r; 0*chnkr.r(1,:,:)];
chnkr.nstor = [chnkr.n; 0*chnkr.r(1,:,:)];
chnkr.dstor = [chnkr.d; 0*chnkr.r(1,:,:)];
chnkr.d2stor = [chnkr.d2; 0*chnkr.r(1,:,:)];
chnkr = sort(chnkr);


figure(1); clf
plot(S,rand(S.npatches,1))
hold on
% plot(chnkr,'x-')
view(0,90)

% Parameters
 

gs_kern_0 = kernel('h','s',1);
gs_kern_0 = gs_kern_0.eval;
% gs_kern_0 = @(s,t) surfwave.flex.kern(s,t,'gs_s',nu,rts_ext, ejs_ext);

gs_kern = gs_kern_0;
% gs_kern = gs_kern_d;
ib2v = 1;

%% Creating geometry and grid


% S = get_surfer(4,8,2*rad);

% ffun = @(t) starfish(t,5,0.,[],[],rad);
% cparams = []; cparams.maxchunklen = 4;
% chnkr = chunkerfunc(ffun,cparams);
% tic;
% S = get_surfer_chnkr(10,8,chnkr);
% tgeom = toc;

ngrid = 1600;
L = 3.5*max(abs(S.r(:)));
xgrid = linspace(-L,L,ngrid+1); 
xgrid = xgrid(1:end-1);
% xgrid = xgrid(2:end);
[xxgrid, yygrid, zzgrid] = meshgrid(xgrid,xgrid,0);

rgrid = [xxgrid(:).'; yygrid(:).'; zzgrid(:).'];

ntarg = 3e4;
rtarg = [2*max(abs(S.r(:)));2*max(abs(S.r(:)));0] .* (rand(3,ntarg) - 0.5)*2;

% thetas = 2*pi*rand(1,ntarg); radtarg = rand(1,ntarg)+15;
% rtarg = radtarg .* [cos(thetas); sin(thetas);0*thetas];

% srcinfo = []; srcinfo.r = S.r; srcinfo.n = [1;0;0] + 0*S.r;
targinfo = []; 
targinfo.r = rtarg; %targinfo.r = S.r;
targinfo.n = randn(3,size(targinfo.r(:,:) ,2)); targinfo.d = randn(3,size(targinfo.r(:,:) ,2)); targinfo.d2 = randn(3,size(targinfo.r(:,:) ,2));
srcinfo = S; 
targinfo = chnkr;
% 
% srcinfo = chnkr; 
targinfo = S;
%% full apply setup

npxy = 9;
npxy = 11;
dx = 2*L/ngrid;
rpxy = 2 * (npxy-1)/2 * (dx);
rbin = 2;


[spread_info,times] = precom_eval_var(gs_kern_0,gs_kern,zk_ext,srcinfo,targinfo,L,dx,ngrid,rgrid,npxy,rpxy,rbin,ib2v,1);


% % [Aspread_s,Aspread_t,Asubtract,kern_hat,rsort,isort,rsort_t,isort_t,times] = precom_eval_var(gs_kern_0,zk_ext,srcinfo,gs_kern_d,targinfo,gs_kern_0,gs_kern_d,L,dx,ngrid,rgrid,npxy,rpxy,rbin,1,1);
% [Aspread_s,Aspread_t,Asubtract,kern_hat,rsort,isort,rsort_t,isort_t,times] = precom_eval_var(gs_kern_0,zk_ext,srcinfo,gs_kern_0,targinfo,gs_kern_d,gs_kern_sp,L,dx,ngrid,rgrid,npxy,rpxy,rbin,0,1);
% % [Aspread_s,Aspread_t,Asubtract,kern_hat,rsort,isort,rsort_t,isort_t,times] = precom_eval_var(gs_kern_0,zk_ext,srcinfo,gs_kern_0,targinfo,gs_kern_0,gs_kern_0,L,dx,ngrid,rgrid,npxy,rpxy,rbin,1,1);


%% full spread apply
inz = randperm(size(srcinfo.r(:,:),2),100);
% inz = 1:size(srcinfo.r(:,:),2);

rsort_t = spread_info.targ_sort;

% inz = 1;
sigma = zeros(size(srcinfo.r(:,:),2),1); sigma(inz) = ones(length(inz),1);
tic;
spread_info.wts = ones(size(srcinfo.r(:,:),2),1);
spread_info.Aquad = sparse(size(targinfo.r(:,:),2),size(srcinfo.r(:,:),2));
u = apply_eval(sigma,spread_info);
% u = apply_eval(sigma,ones(S.npts,1),Aspread_s,Aspread_t,Asubtract,kern_hat,isort,isort_t,ngrid, sparse(size(rtarg,2),S.npts));
u = u(spread_info.isort_t);
toc;

tic;
uplot = gs_kern(struct('r',srcinfo.r(:,inz),'n',srcinfo.n(:,inz)),rsort_t)*sigma(inz);

% uplot = gs_kern_0(struct('r',srcinfo.r(:,inz)),struct('r',rsort_t))*sigma(inz);
% uplot = gs_kern_d(struct('r',srcinfo.r(:,inz),'n',srcinfo.n(:,inz)),rsort_t)*sigma(inz);
% uplot = gs_kern_sp(struct('r',srcinfo.r(:,inz),'n',srcinfo.n(:,inz)),rsort_t)*sigma(inz);
% u = u*norm(uplot,inf)/norm(u,inf);
toc
figure(5)
subplot(1,3,1)
my_scatter(rsort_t.r,[],real(u/norm(uplot,inf)),'.'); colorbar
hold on
my_scatter(srcinfo.r(:,inz),'x','linewidth',2)
hold off

disp('123')

% uplot = gs_kern(struct('r',S.r(:,1)),struct('r',rsort));
subplot(1,3,2)
my_scatter(rsort_t.r,[],real(uplot/norm(uplot,inf)),'.');colorbar
subplot(1,3,3)
my_scatter(rsort_t.r,[],log10(abs(u-uplot)/norm(uplot,inf)),'.'); colorbar
% clim([-12,-7])
% hold on, plot(chnkr,'linewidth',2), hold off

tic;
str = spread_info.Aspread_s*sigma(spread_info.isort); str = full(str);
tspread = toc

tic;
str_hat = fft2(reshape(str,ngrid,ngrid));
u_hat = spread_info.kern_hat.*str_hat;
ugrid = ifft2(u_hat);
tconv = toc

tic;
u2 = spread_info.Aspread_t.'*ugrid(:);
tinterp = toc
%
figure(4)
subplot(1,2,1)
my_scatter(rsort_t.r,[],log10(abs(u-uplot)/norm(uplot,inf)),'o','linewidth',2); colorbar
% clim([-12,-7])

subplot(1,2,2)
my_scatter(rsort_t.r,[],log10(abs(u2-uplot)/norm(uplot,inf)),'.'); colorbar
% clim([-12,-7])




%%
% sigma = exp(-vecnorm(S.r).^2/100).';
% t1 = tic;
% u = apply_eval(sigma,S.wts,Aspread_s,Aspread_t,Asubtract,kern_hat,isort,isort_t,ngrid, sparse(size(rtarg,2),S.npts));
% toc(t1);
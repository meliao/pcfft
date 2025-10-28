addpath(genpath('~/FastAlgs/chunkie'))
addpath(genpath('~/FastAlgs2/fmm3dbie_log_quad'))
rng(123)

%% set parameters
alpha_int = 4;
beta_int = 1;
[rts_int,ejs_int] = surfwave.capillary.find_roots(alpha_int,beta_int,1);

% alpha_ext = 10;
% beta_ext = 1;
alpha_ext = 1;
beta_ext = 1200;
[rts_ext,ejs_ext] = surfwave.capillary.find_roots(alpha_ext,beta_ext,1);

zpars = [rts_ext;ejs_ext];

zk_int = rts_int(1);
zk_ext = rts_ext(1);

gs_kern = @(s,t) surfwave.kern(rts_ext,ejs_ext,s,t,'gs_s');
% gs_kern = @(s,t) surfwave.kern(rts_ext,ejs_ext,s,t,'gphi_s');

%% Creating geometry and grid


% rad = 1;
% S = geometries.disk([rad,rad],[],[4 3 6],12);


load +geometries/srfmat.mat
S = rotate(srfr,[pi 0 0]);
S = translate(S,-mean(S.r,2));
% S = oversample(S,12);

% ffun = @(t) starfish(t,5,0.3,[],[],rad);
% cparams = []; cparams.maxchunklen = 4;
% chnkr = chunkerfunc(ffun,cparams);
% tic;
% S = get_surfer_chnkr(10,8,chnkr);
% tgeom = toc;

% return

ngrid = 2*800;
L = 2.5*max(abs(S.r(:)));
xgrid = linspace(-L,L,ngrid+1); 
xgrid = xgrid(1:end-1);
% xgrid = xgrid(2:end);
[xxgrid, yygrid, zzgrid] = meshgrid(xgrid,xgrid,0);

rgrid = [xxgrid(:).'; yygrid(:).'; zzgrid(:).'];

%% full apply setup

npxy = 9;
% npxy = 11;
dx = 2*L/ngrid;
rpxy = 2 * (npxy-1)/2 * (dx);
rbin = 2;
eps = 1e-10;


[Aspread,Asubtract,kern_hat,rsort,isort] = precom_apply(gs_kern,zk_ext,S,L,dx,ngrid,rgrid,npxy,rpxy,rbin,1);
% tic; Aquad =surfwave.capillary.get_quad_cor(S,eps,zpars,0);toc;
% 
% tic;
% [i,j,vals] = find(Aquad);
% rsrc = S.r(:,i);
% rtarg = S.r(:,j);
% rnear = rtarg - rsrc;
% vals_sub = gs_kern(struct('r',[0;0;0]),struct('r',rnear));
% Aquad_sub = sparse(i,j, vals-vals_sub);
% toc;
        



%% full spread apply
inz = randperm(S.npts,10);
% inz = 1;
sigma = zeros(S.npts,1); sigma(inz) = ones(length(inz),1);

tic
u = apply_mat(sigma./S.wts,S.wts,Aspread,Asubtract,kern_hat,isort,ngrid,0);
% u(isort) = u;
u = u(isort);
% u = u +Aquad_sub*sigma;
tapply = toc

uplot = gs_kern(struct('r',S.r(:,inz)),struct('r',rsort))*sigma(inz);
figure(5)
subplot(1,3,1)
my_scatter(rsort,[],real(u/norm(uplot,inf)),'.'); colorbar
hold on
my_scatter(S.r(:,inz),'x','linewidth',2)
hold off

% uplot = gs_kern(struct('r',S.r(:,1)),struct('r',rsort));
subplot(1,3,2)
my_scatter(rsort,[],real(uplot/norm(uplot,inf)),'.');colorbar
subplot(1,3,3)
my_scatter(rsort,[],log10(abs(u-uplot)/norm(uplot,inf)),'.'); colorbar
clim([-12,-7])
% hold on, plot(chnkr,'linewidth',2), hold off

tic;
str = Aspread*sigma(isort); str = full(str);
% tspread = toc

tic;
str_hat = fft2(reshape(str,ngrid,ngrid));
u_hat = kern_hat.*str_hat;
ugrid = ifft2(u_hat);
% tconv = toc

tic;
u2 = Aspread.'*ugrid(:);

u2 =u2 + Asubtract*sigma(isort);
% tinterp = toc
figure(4)
subplot(1,2,1)
my_scatter(rsort,[],log10(abs(u-uplot)/norm(uplot,inf)),'.'); colorbar
% clim([-12,-7])

subplot(1,2,2)
my_scatter(rsort,30,log10(abs(u2-uplot)/norm(uplot,inf)),'.'); colorbar
% my_scatter(rsort,[],log10(abs(u2)/norm(uplot,inf)),'.'); colorbar
% clim([-12,-7])

max(abs(u-uplot))


% %%
% sigma = exp(-vecnorm(S.r).^2/100).';
% t1 = tic;
% u = apply_mat(sigma,S.wts,Aspread,Asubtract,kern_hat,isort,ngrid,Aquad_sub);
% toc(t1);
% 
% Amat = @(x) x + apply_mat(x,S.wts,Aspread,Asubtract,kern_hat,isort,ngrid,Aquad_sub);
% 
% tic;
% x = gmres(Amat,sigma,[],1e-6,1000);
% toc;




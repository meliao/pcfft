
rad = pi/2;
beta = 0.5;
gamma = 1;
[rts,ejs] = surfwave.capillary.find_roots(beta,gamma);

zpars = [rts;ejs];

zk_ext = rts(1)

% Creating geometry

S = geometries.disk([rad,rad],[],[3 3 3],6);
% S = geometries.disk([rad,rad],[],[5 5 5],6);

maxrad = max(S.r(:));

norder = S.norders(1);

% Creating geometry

nch = 4*3;
cparams.ta = pi/nch;
cparams.tb = 2*pi + cparams.ta;
chnkr = chunkerfuncuni(@(t) ellipse(t),nch,cparams);
chnkr = chnkr.move([0;0],[0;0],0,rad);
chnkr = sort(chnkr);


gamma = 0.6;
ngrid = 2*1600; % should scale with density of pts
L = 2.5*max(abs(S.r(:)));
xgrid = linspace(-L,L,ngrid+1); 
xgrid = xgrid(1:end-1);
% xgrid = xgrid(2:end);
[xxgrid, yygrid, zzgrid] = meshgrid(xgrid,xgrid,0);

rgrid = [xxgrid(:).'; yygrid(:).'; zzgrid(:).'];

npxy = 11; % no le toque
dx = 2*L/ngrid;
rpxy = 2 * (npxy-1)/2 * (dx);
% rpxy = 3 * (npxy-1)/2 * (dx);
rbin = 2; % no le toque

grid_param.L = L;
grid_param.dx = dx;
grid_param.ngrid = ngrid;
grid_param.rgrid = rgrid;
grid_param.npxy = npxy;
grid_param.rpxy = rpxy;
grid_param.rbin = rbin;


eps = 1e-9;

gs_kern = @(s,t) surfwave.capillary.kern(rts,ejs,s,t,'gs_s');
gp_kern = @(s,t) surfwave.capillary.kern(rts,ejs,s,t,'gphi_s');
% lap_gp_kern = @(s,t) surfwave.capillary.kern(rts,ejs,s,t,'lap_gphi');
s3d_gp_kern = @(s,t) surfwave.capillary.kern(rts,ejs,s,t,'s3d_gphi');
s3d_kern = @(s,t) surfwave.flex.lap3dkern(s.r,t.r);
s3d_plus_s3d_gp_kern = @(s,t) s3d_kern(s,t) + s3d_gp_kern(s,t);
lap_gp_kern = @(s,t) -2/beta*s3d_kern(s,t) ... 
    - 2/beta*s3d_gp_kern(s,t) + gamma/beta*gp_kern(s,t);

% tic;
% spread_info = precom_eval_var(gs_kern,gs_kern,zk_ext,S,S,L,dx,ngrid,rgrid,npxy,rpxy,rbin,1,0);
% toc;
% t1 = tic; 
% Gsquad_sub = surfwave.capillary.get_quad_cor_sub(S,gs_kern,eps,zpars,0);
% tquad = toc(t1)
% 
% %%
% 
% spread_info.wts = S.wts;
% spread_info.Aquad = Gsquad_sub;
% 
% 
% Gspread_S = spread_info.Aspread_s;
% Gskern_hat = spread_info.kern_hat;
% Srsort = spread_info.src_sort;
% Sisort = spread_info.isort;
% Sisortinv = 1:S.npts;
% Sisortinv(Sisort) = Sisortinv;
% 
% Gs_addsub_S = spread_info.Asubtract(Sisortinv,Sisortinv).*S.wts(:).' + Gsquad_sub;
% 
% sigma = randn(S.npts,1);

%%

% 
% t1 = tic;
% spread_info = precom_eval_var(gs_kern,gs_kern,zk_ext,S,S,L,dx,ngrid,rgrid,npxy,rpxy,rbin,1,0);
% toc(t1)
% t1 = tic; 
% Gsquad_sub = surfwave.capillary.get_quad_cor_sub(S,gs_kern,eps,zpars,0);
% tquad = toc(t1)
% 
% Gsspread_S = spread_info.Aspread_s;
% Gskern_hat = spread_info.kern_hat;
% Srsorted = spread_info.src_sort;
% Sisort = spread_info.isort;
% Sisortinv = 1:S.npts;
% Sisortinv(Sisort) = Sisortinv;
% Gs_addsub_S = spread_info.Asubtract(Sisortinv,Sisortinv).*S.wts(:).' + Gsquad_sub;
% 
% 
% spread_info.wts = S.wts; spread_info.Aquad = Gsquad_sub;
% gs_apply = @(mu) apply_eval(mu,spread_info);
% 
% 
% t1 = tic;
% spread_info = precom_eval_var(gp_kern,gp_kern,zk_ext,S,S,L,dx,ngrid,rgrid,npxy,rpxy,rbin,1,0);
% toc(t1)
% t1 = tic; 
% Gpquad_sub = surfwave.capillary.get_quad_cor_sub(S,gp_kern,eps,zpars,1);
% tquad = toc(t1)
% 
% Gpspread_S = spread_info.Aspread_s;
% Gpkern_hat = spread_info.kern_hat;
% Gp_addsub_S = spread_info.Asubtract(Sisortinv,Sisortinv).*S.wts(:).' + Gpquad_sub;
% 
% spread_info.wts = S.wts; spread_info.Aquad = Gpquad_sub;
% gp_apply = @(mu) apply_eval(mu,spread_info);
% 
% % S3d plus S3d Gphi
% 
% t1 = tic;
% spread_info = precom_eval_var(s3d_plus_s3d_gp_kern,s3d_plus_s3d_gp_kern,zk_ext,S,S,L,dx,ngrid,rgrid,npxy,rpxy,rbin,1,0);
% toc(t1);
% t1 = tic; 
% s3d_quad_sub = surfwave.capillary.get_quad_cor_sub(S,s3d_kern,eps,zpars,6);
% s3d_gphi_quad_sub = surfwave.capillary.get_quad_cor_sub(S,s3d_gp_kern,eps,zpars,5);
% tquad = toc(t1)
% 
% Gs3d_plus_s3d_gpspread_S = spread_info.Aspread_s;
% Gs3d_plus_s3d_gpkern_hat = spread_info.kern_hat;
% Gs3d_plus_s3d_gp_addsub_S = spread_info.Asubtract(Sisortinv,Sisortinv).*S.wts(:).' + s3d_quad_sub + s3d_gphi_quad_sub;
% 
% spread_info.wts = S.wts; spread_info.Aquad = s3d_quad_sub + s3d_gphi_quad_sub;
% s3d_plus_s3d_gp_apply = @(mu) apply_eval(mu,spread_info);
% 
% lap_gp_apply = @(mu) -2/beta*s3d_plus_s3d_gp_apply(mu) + gamma/beta*gp_apply(mu);
% 
% tv2v = toc(t4)
% 
% 
% 
% 
% 
% 
% %%
% sigma = randn(S.npts,1);
% 
% 
% sigmasort = sigma(spread_info.isort).*spread_info.wts(spread_info.isort);
% str = Gsspread_S*sigmasort; str = full(str);
% 
% str_hat = fft2(reshape(str,spread_info.ngrid,spread_info.ngrid));
% u_hat = Gskern_hat.*str_hat;
% ugrid = ifft2(u_hat);
% 
% u = Gs_addsub_S*sigma +  Gsspread_S(:,Sisortinv).'*ugrid(:);
% 
% norm(gs_apply(sigma) - u)
% 
% 
% 
% str = Gpspread_S*sigmasort; str = full(str);
% 
% str_hat = fft2(reshape(str,spread_info.ngrid,spread_info.ngrid));
% u_hat = Gpkern_hat.*str_hat;
% ugrid = ifft2(u_hat);
% 
% u = Gp_addsub_S*sigma +  Gpspread_S(:,Sisortinv).'*ugrid(:);
% 
% norm(gp_apply(sigma) - u)
% 
% 
% 
% str = Gs3d_plus_s3d_gpspread_S*sigmasort; str = full(str);
% 
% str_hat = fft2(reshape(str,spread_info.ngrid,spread_info.ngrid));
% u_hat = Gs3d_plus_s3d_gpkern_hat.*str_hat;
% ugrid = ifft2(u_hat);
% 
% u2 = Gs3d_plus_s3d_gp_addsub_S*sigma +  Gs3d_plus_s3d_gpspread_S(:,Sisortinv).'*ugrid(:);
% 
% norm(lap_gp_apply(sigma) - (-2/beta*u2 + gamma/beta*u))
% 

%%
gs_kern = @(s,t) surfwave.capillary.kern(rts,ejs,s,t,'gs_s');
gp_kern = @(s,t) surfwave.capillary.kern(rts,ejs,s,t,'gphi_s');
% lap_gp_kern = @(s,t) surfwave.capillary.kern(rts,ejs,s,t,'lap_gphi');
s3d_gp_kern = @(s,t) surfwave.capillary.kern(rts,ejs,s,t,'s3d_gphi');
s3d_kern = @(s,t) surfwave.flex.lap3dkern(s.r,t.r);
s3d_plus_s3d_gp_kern = @(s,t) s3d_kern(s,t) + s3d_gp_kern(s,t);
lap_gp_kern = @(s,t) -2/beta*s3d_kern(s,t) ... 
    - 2/beta*s3d_gp_kern(s,t) + gamma/beta*gp_kern(s,t);

%Gs
spread_info = precom_eval_var(gs_kern,gs_kern,zk_ext,S,S,L,dx,ngrid,rgrid,npxy,rpxy,rbin,1,0);
t1 = tic; 
Gsquad_sub = surfwave.capillary.get_quad_cor_sub(S,gs_kern,eps,zpars,0);
tquad = toc(t1)

spread_info.wts = S.wts; spread_info.Aquad = Gsquad_sub;
gs_apply = @(mu) apply_eval(mu,spread_info);


%%


spread_info = precom_eval_var(gp_kern,gp_kern,zk_ext,S,S,L,dx,ngrid,rgrid,npxy,rpxy,rbin,1,0);
t1 = tic; 
Gpquad_sub = surfwave.capillary.get_quad_cor_sub(S,gp_kern,eps,zpars,1);
tquad = toc(t1)

spread_info.wts = S.wts; spread_info.Aquad = Gpquad_sub;
gp_apply = @(mu) apply_eval(mu,spread_info);

% S3d plus S3d Gphi

spread_info = precom_eval_var(s3d_plus_s3d_gp_kern,s3d_plus_s3d_gp_kern,zk_ext,S,S,L,dx,ngrid,rgrid,npxy,rpxy,rbin,1,0);
t1 = tic; 
s3d_quad_sub = surfwave.capillary.get_quad_cor_sub(S,s3d_kern,eps,zpars,6);
s3d_gphi_quad_sub = surfwave.capillary.get_quad_cor_sub(S,s3d_gp_kern,eps,zpars,5);
tquad = toc(t1)

spread_info.wts = S.wts; spread_info.Aquad = s3d_quad_sub + s3d_gphi_quad_sub;
s3d_plus_s3d_gp_apply = @(mu) apply_eval(mu,spread_info);

lap_gp_apply = @(mu) -2/beta*s3d_plus_s3d_gp_apply(mu) + gamma/beta*gp_apply(mu);


%% Getting for S_S and S_\phi (boundary -> volume)

opts = [];
opts.eps = eps;
opts.corrections = true;

tic
gs_b2v_corr = chunkerkernevalmat(chnkr,gs_kern,S.r(1:2,:),opts);
gp_b2v_corr = chunkerkernevalmat(chnkr,gp_kern,S.r(1:2,:),opts);
t5 = toc

srcinfo = chnkr2; 
targinfo = S;

[spread_info,times] = precom_eval_var(gs_kern,gs_kern,zk_ext,srcinfo,targinfo,L,dx,ngrid,rgrid,npxy,rpxy,rbin,1,0);
spread_info.wts = chnkr2.wts(:); spread_info.Aquad = gs_b2v_corr;
gs_b2v_apply = @(rho) apply_eval(rho,spread_info);

[spread_info,times] = precom_eval_var(gp_kern,gp_kern,zk_ext,srcinfo,targinfo,L,dx,ngrid,rgrid,npxy,rpxy,rbin,1,0);
spread_info.wts = chnkr2.wts(:); spread_info.Aquad = gp_b2v_corr;
gp_b2v_apply = @(rho) apply_eval(rho,spread_info);


%% Forming matrices for V' (volume -> boundary)

gs_sp_kern = @(s,t) surfwave.capillary.kern(rts, ejs, s, t,'gs_sprime');
gp_sp_kern = @(s,t) surfwave.capillary.kern(rts, ejs, s, t,'gphi_sprime');

% Forming matrix for V' (volume -> boundary)

targinfo = zeros(11,chnkr.npt);
targinfo(1:2,:) = chnkr.r(:,:);
targinfo(10:11,:) = chnkr.n(:,:);

% eps = 1e-9;

% tic
% gsgradv2b = surfwave.capillary.v2bmat(S,targinfo,zpars, eps);
% t3 = toc
% 
% gsgradv2b = reshape(gsgradv2b, [S.npts chnkr.npt]).';

ekerns(2,1) = kernel();
ekerns(1) = gs_sp_kern;
ekerns(2) = gp_sp_kern;

tic
[gs_v2b_corr,gp_v2b_corr] = surfwave.capillary.get_quad_cor_v2b(S,targinfo,ekerns,eps,zpars); %,uv_bndry,patch_id);
t3 = toc

nover = S.norders(1) + 4;
srfr_over = oversample(S,nover);

srcinfo = srfr_over;
targinfo = chnkr2;

[spread_info,times] = precom_eval_var(gs_kern,gs_sp_kern,zk_ext,srcinfo,targinfo,L,dx,ngrid,rgrid,npxy,rpxy,rbin,0,0);
spread_info.wts = srfr_over.wts;
gs_sp_v2b_apply = @(mu) apply_eval_over(S,mu,spread_info,gs_v2b_corr);

[spread_info,times] = precom_eval_var(gp_kern,gp_sp_kern,zk_ext,srcinfo,targinfo,L,dx,ngrid,rgrid,npxy,rpxy,rbin,0,0);
spread_info.wts = srfr_over.wts;
spread_info.Aquad = gp_v2b_corr;
gp_sp_v2b_apply = @(mu) apply_eval_over(S,mu,spread_info,gp_v2b_corr);
% gp_sp_v2b_apply = @(mu) apply_eval(mu,spread_info);

%% Forming matrix for S prime and D prime (boundary -> boundary)

tic
Sprimeb2b = chunkermat(chnkr,gs_sp_kern);
t6 = toc

%%


b2v_apply = @(rho) 1/2*gs_b2v_apply(rho) - gp_b2v_apply(rho);
v2v_apply = @(mu) 1/2*mu + 1/4*(abs(gamma)-gamma)*gs_apply(mu) + (1-abs(gamma))/2*gp_apply(mu) + beta/2*lap_gp_apply(mu);
b2b_apply =  -1/(2*beta)*eye(chnkr.npt) + 1/2*Sprimeb2b;
v2b_apply = @(mu) 1/4*(abs(gamma)-gamma)*gs_sp_v2b_apply(mu)+1/2*gp_sp_v2b_apply(mu);

%% use new Capillary wrapper
norderup = 4;
 [Gsspread_S,Gsspread_Sover,Gskern_hat,Gs_addsub_S,...
    Gpspread_S,Gpspread_Sover, Gpkern_hat, Gp_addsub_S,...
    Gs3d_plus_s3d_gpspread_S,Gs3d_plus_s3d_gpspread_Sover,Gs3d_plus_s3d_gpkern_hat,Gs3d_plus_s3d_gp_addsub_S,...
    Gsspread_C,Gs_addsub_C2S,Gpspread_C,Gp_addsub_C2S,...
    Gsprimespread_C,Gsprime_addsub_S2C,Gpprimespread_C,Gpprime_addsub_S2C,...
    Sprimeb2b,val2over,times] = build_capillary_neu(S,chnkr,rts,ejs,beta,gamma,zpars,zk_ext,eps,grid_param,norderup,[],[]);


%% use new Capillary apply wrapper
nover = S.norders(1) + norderup; Sover = oversample(S,nover);
Swts = Sover.wts(:); Snpt = S.npts; Cwts = chnkr.wts(:); 

A_app = @(dense) apply_capillary_neu(beta,gamma,Gsspread_S,Gsspread_Sover,Gskern_hat,Gs_addsub_S,...
    Gpspread_S,Gpspread_Sover, Gpkern_hat, Gp_addsub_S,...
    Gs3d_plus_s3d_gpspread_S,Gs3d_plus_s3d_gpspread_Sover,Gs3d_plus_s3d_gpkern_hat,Gs3d_plus_s3d_gp_addsub_S,...
    Gsspread_C,Gs_addsub_C2S,Gpspread_C,Gp_addsub_C2S,...
    Gsprimespread_C,Gsprime_addsub_S2C,Gpprimespread_C,Gpprime_addsub_S2C,...
    Sprimeb2b,val2over,Swts,Cwts,dense,Snpt);

% dense = [randn(S.npts,1);randn(chnkr.npt,1)];
% dense = [ones(S.npts,1);ones(chnkr.npt,1)];
% dense = [exp(S.r(1,:)).';exp(chnkr.r(1,:)).'];
dense = [ones(S.npts,1);zeros(chnkr.npt,1)];



mu = dense(1:S.npts);
rho = dense(Snpt+1:end);

% u = gs_apply(mu) + gp_apply(mu) + lap_gp_apply(mu);
% u =gp_apply(mu);
% u = v2v_apply(mu);
% b = v2b_apply(mu);
% 
% u = b2v_apply(rho);
% b = b2b_apply*rho;

u = v2v_apply(mu) + b2v_apply(rho);
b = v2b_apply(mu) + b2b_apply*rho;

uout = A_app(dense);
u2 = uout(1:S.npts);
b2 = uout(Snpt+1:end);

norm(u-u2)
norm(b-b2)




%% set parameters
% alpha_int = 4;
% beta_int = 1;
% [rts_int,ejs_int] = surfwave.capillary.find_roots(alpha_int,beta_int,1);
% 
% % alpha_ext = 4;
% % beta_ext = 1;
% alpha_ext = 1;
% beta_ext = 1200;
% [rts_ext,ejs_ext] = surfwave.capillary.find_roots(alpha_ext,beta_ext,1);

% alpha_ext = 1.5;
beta_ext = 1.5;
gamma_ext = 1;
nu = 0.3;
% [rts_ext,ejs_ext] = surfwave.flex.find_roots(alpha_ext,beta_ext,gamma_ext);
[rts_ext,ejs_ext] = surfwave.flex.find_roots(beta_ext,gamma_ext);

zpars = [rts_ext;ejs_ext];

% zk_int = rts_int(1);
zk_ext = rts_ext(1);

% gs_kern = @(s,t) surfwave.flex.kern(rts_ext,ejs_ext,s,t,'gs_s');
gs_kern = kernel(@(s,t) surfwave.flex.kern(s,t,'gs_s',nu,rts_ext,ejs_ext));
% gs_kern = @(s,t) surfwave.flex.kern(s,t,'free plate eval second gs',nu,rts_ext,ejs_ext);
% gs_kern = @(s,t) surfwave.flex.kern(s,t,'free plate eval second gphi',nu,rts_ext,ejs_ext);

%% Creating geometry and grid


% rad = 30; rad = 2*pi/zk_ext*10;
% S = geometries.disk([rad,rad],[],[4 3 4],9);
% 
% ngrid = 100;
% L = 1.5*max(abs(S.r(:)));
% xgrid = linspace(-L,L,ngrid+1); xgrid = xgrid(1:end-1);
% [xxgrid, yygrid, zzgrid] = meshgrid(xgrid,xgrid,0);
% 
% rgrid = [xxgrid(:).'; yygrid(:).'; zzgrid(:).'];


% %% setup proxy
% 
% 
% npxy = 9;
% 
% % dxgrid = 2*L/ngrid;
% dxgrid = 2*pi/zk_ext/100;
% rpxy = 1.5 * (npxy-1)/2 * (dxgrid);
% % rpxy = 2 * (npxy-1)/2 * (dxgrid);
% [pxy,str2pxy,pxygrid] = setup_proxy(dxgrid,npxy,rpxy,gs_kern,1e-10,zk_ext);
% % [pxy,str2pxy,pxygrid] = setup_proxy(dxgrid,npxy,rpxy,gs_kern);
% rng(123)
% rpxy = 1.5*rpxy;
% 
% figure(1);clf
% % plot(S), hold on
% my_scatter(pxy,'x')
% hold on
% my_scatter(pxygrid,'.')
% % my_scatter(rgrid,'.')
% % my_scatter(S.r,'x')
% hold off
% 
% 
% src = [dxgrid/3;dxgrid/2;0];
% src2pxy = gs_kern(struct('r', src), struct('r',pxy));
% 
% str = str2pxy\src2pxy;
% 
% theta = 0.1;
% targ = rpxy*1.3 * [cos(theta); sin(theta); 0];
% % targ = rpxy*1.1 * [cos(theta); sin(theta); 0];
% 
% src2targ = gs_kern(struct('r', src), struct('r',targ));
% grid2targ = gs_kern(struct('r', pxygrid), struct('r',targ));
% 
% hold on, my_scatter(targ,'x', 'linewidth',2), hold off
% 
% norm(src2targ - grid2targ*str)
% 
% ntarg = 2e4;
% rs = 0.2*rpxy*rand(1,ntarg)+rpxy;
% thetas = 2*pi*rand(1,ntarg);
% targs = rs.* [cos(thetas); sin(thetas); 0*thetas];
% 
% nsrc = 100;
% src = 2*[dxgrid;dxgrid;0].*(rand(3,nsrc)-0.5);
% 
% tic;
% src2pxy = gs_kern(struct('r', src), struct('r',pxy));
% src2targ = gs_kern(struct('r', src), struct('r',targs));
% grid2targ = gs_kern(struct('r', pxygrid), struct('r',targs));
% toc;
% 
% 
% errs = vecnorm(src2targ- grid2targ*(str2pxy\src2pxy),Inf,2);
% % %%
% figure(2);clf
% my_scatter(targs,[],log10(errs),'.')
% hold on
% my_scatter(src,'.')
% my_scatter(pxygrid,'o')
% my_scatter(pxy,'r^','linewidth',2)
% hold off
% colorbar
% title('log10 max proxy error')
% % return
% %%
% % rpxy = 2 * (npxy-1)/2 * (dxgrid);
% % [pxy,str2pxy,pxygrid] = setup_proxy(dxgrid,npxy,rpxy,gs_kern);
% ntarg = 0.5e4;
% % rs = 2*rpxy*rand(1,ntarg)+rpxy;
% rs = 2*(2*pi/zk_ext)*rand(1,ntarg)+rpxy;
% thetas = 2*pi*rand(1,ntarg);
% targs = rs.* [cos(thetas); sin(thetas); 0*thetas];
% 
% nsrc = 100;
% src = 3*[dxgrid;dxgrid;0].*(rand(3,nsrc)-0.5);
% 
% tic;
% src2pxy = gs_kern(struct('r', src), struct('r',pxy));
% src2targ = gs_kern(struct('r', src), struct('r',targs));
% grid2targ = gs_kern(struct('r', pxygrid), struct('r',targs));
% toc;
% 
% 
% errs = vecnorm(src2targ- grid2targ*(str2pxy\src2pxy),Inf,2);
% figure(3);clf; t= tiledlayout('flow'); t.TileSpacing = 'tight';
% nexttile()
% my_scatter(targs,[],log10(errs),'.')
% hold on
% my_scatter(src,'.')
% my_scatter(pxygrid,'o')
% % my_scatter(pxy,'r^','linewidth',2)
% hold off
% colorbar
% title('log10 max proxy error')
% clim([-7,-6])
% 
% 
% nexttile()
% my_scatter(targs,[],log10(vecnorm(src2targ,Inf,2)),'.')
% hold on
% my_scatter(src,'.')
% my_scatter(pxygrid,'o')
% % my_scatter(pxy,'r^','linewidth',2)
% hold off
% colorbar
% title('log10 max kern')
% 
% 
% max(errs)


%%
npxy = 13;


% dxs = 2*pi./zk_ext./[2,4,10.^(1:3)];
dxs = 2*pi./zk_ext./[4,10.^(1:3)];
% tols = 10.^(-7:-1:-12);
tols = 10.^(-11:-1:-13);
% tols = 1e-12;


errsmax = zeros(length(dxs), length(tols));
nprox = zeros(length(dxs), length(tols));



for i = 1:length(dxs)
    tic;
    dxgrid = dxs(i);
    rpxy = 3 * (npxy-1)/2 * (dxgrid);

   % rpxy = 1.*rpxy;
    ntarg = 0.5e4;
    % rs = 2*rpxy*rand(1,ntarg)+rpxy;
    rs = 4*(2*pi/zk_ext)*rand(1,ntarg)+rpxy;
    % rs = 1*(2*pi/zk_ext)*rand(1,ntarg)+rpxy;
    thetas = 2*pi*rand(1,ntarg);
    targs = rs.* [cos(thetas); sin(thetas); 0*thetas];

    nsrc = 100;
    src = 2*[dxgrid;dxgrid;0].*(rand(3,nsrc)-0.5);
    src2targ = gs_kern.eval(struct('r', src), struct('r',targs));

    
    for j = 1:length(tols)
    
    [pxy,str2pxy,pxygrid] = setup_proxy(dxgrid,npxy,rpxy,gs_kern,tols(j));
    [pxy,str2pxy,pxygrid] = setup_proxy(dxgrid,npxy,rpxy,gs_kern,tols(j),zk_ext);
    
    src2pxy = gs_kern.eval(struct('r', src), struct('r',pxy));
    
    grid2targ = gs_kern.eval(struct('r', pxygrid), struct('r',targs));
    
    % errs = vecnorm(src2targ- grid2targ*(str2pxy*src2pxy),Inf,2);
    errs = vecnorm(src2targ- grid2targ*(str2pxy\src2pxy),Inf,2);
    errsmax(i,j) = max(errs);
    nprox(i,j) = size(pxy,2);
    end
    toc;
end

errsmax
nprox
% %%
figure(7);clf
plot(log10(dxs), log10(errsmax),'o-','LineWidth',2)
xlabel('$\log_{10} dx$','interpreter','latex')
ylabel('$\log_{10}$ proxy error','interpreter','latex')
legend(num2str(tols'),'interpreter','latex')
set(gca,'fontsize',14)

% %%
figure(1);clf
% plot(S), hold on
my_scatter(pxy,'x')
hold on
my_scatter(pxygrid,'.')
my_scatter(src,'o')
% my_scatter(targs,'.')
% my_scatter(S.r,'x')
hold off
axis equal


% %%
figure(3);clf
my_scatter(targs,[],log10(errs).','.')
hold on
my_scatter(pxy,'x')

colorbar







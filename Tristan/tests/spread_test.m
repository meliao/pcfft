
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


rad = 1;
S = geometries.disk([rad,rad],[],2*[4 3 6],8);

ngrid = 200;
L = 2.5*max(abs(S.r(:)));
xgrid = linspace(-L,L,ngrid+1); 
xgrid = xgrid(1:end-1);
% xgrid = xgrid(2:end);
[xxgrid, yygrid, zzgrid] = meshgrid(xgrid,xgrid,0);

rgrid = [xxgrid(:).'; yygrid(:).'; zzgrid(:).'];


% %% setup proxy
% 
% npxy = 9;
dx = 2*L/ngrid;
% rpxy = 2 * (npxy-1)/2 * (dx);
% [pxy,str2pxy,pxygrid] = setup_proxy(dx,npxy,rpxy,gs_kern,1e-8,zk_ext);
% 
% 
rbin = 3;
nbin = floor(ngrid/rbin);
% 
% 
% sigma = randn(S.npts,1);
sigma = zeros(S.npts,1); sigma(1) = 1;
% sigma = exp(-vecnorm(S.r).^2/100).';
% 
% 
% %% test binning
[idstart,isort,rsort] = bin_pts(S.r,rbin,dx,nbin,L);

% str = zeros(ngrid*ngrid,1);

% truebin = floor(1*nbin/6):ceil(5*nbin/6);
% bin_ctrs = zeros(3,nbin*nbin);
% for i = 1:nbin
%     xcenter = rbin * (i-1) * dx - L;
%     for j = 1:nbin
%         ycenter =  rbin * (j-1) * dx - L;
%         rcenter = [xcenter;ycenter;0];
%         bin_ctrs(:,(i-1)*nbin + j)= rcenter;
%     end
% end
% % 
% figure(1);clf
% my_scatter(bin_ctrs,'o')
% hold on
% my_scatter(rgrid,'.')
% my_scatter(S.r,'x')
% hold off
% 
% 
% figure(2); clf
% ibin = 200;
% my_scatter(rsort(:,idstart(ibin):(idstart(ibin+1)-1)))
% hold on
% my_scatter(bin_ctrs(:,ibin),'x')
% hold off
% 
% 
% % verify that no point is too far from its bin center
% rad = 0;
% for ibin = 1:nbin*nbin
%     if idstart(ibin) == idstart(ibin+1), continue, end
%     tmp = max(vecnorm(rsort(:,idstart(ibin):(idstart(ibin+1)-1)) - bin_ctrs(:,ibin),inf));
%     rad = max(tmp,rad);
% end
% 
% [rad,rbin*dx/2,rbin*dx/2-rad]
% 
% %% spread
% tic;
% sigma_sort = sigma(isort);
% 
% str = zeros(ngrid*ngrid,1);
% 
% for ibin = 1:nbin*nbin
%     if idstart(ibin) == idstart(ibin+1), continue, end
% 
%     id_loc = idstart(ibin):(idstart(ibin+1)-1);
%     r_loc = rsort(:,id_loc) - bin_ctrs(:,ibin);
%     sig_loc = sigma_sort(id_loc);
% 
%     loc2pxy = gs_kern(struct('r',r_loc),struct('r',pxy))*sig_loc;
% 
%     str_loc = str2pxy\loc2pxy;
% 
%     id_bin_y = mod(ibin,nbin)+1;
%     id_bin_x = ceil(ibin/nbin);
% 
%     id_grid_x = (id_bin_x*rbin-1) + (-((npxy-1)/2):((npxy-1)/2));
%     id_grid_y = (id_bin_y*rbin-1) + (-((npxy-1)/2):((npxy-1)/2));
%     id_grid_loc = (id_grid_x(:)-1)*ngrid + id_grid_y(:).';
% 
%     str(id_grid_loc(:)) = str(id_grid_loc(:)) + str_loc;
% end
% 
% toc;
% %%
% tic;
% r_loc = 0*rsort;
% for ibin = 1:nbin*nbin
%     if idstart(ibin) == idstart(ibin+1), continue, end
% 
%     id_loc = idstart(ibin):(idstart(ibin+1)-1);
%     r_loc(:,id_loc) = rsort(:,id_loc) - bin_ctrs(:,ibin);
% end
% 
% pts2pxy = gs_kern(struct('r',r_loc),struct('r',pxy));
% str2 = zeros(ngrid*ngrid,1);
% 
% for ibin = 1:nbin*nbin
%     if idstart(ibin) == idstart(ibin+1), continue, end
% 
%     id_loc = idstart(ibin):(idstart(ibin+1)-1);
%     sig_loc = sigma_sort(id_loc);
% 
%     str_loc = str2pxy\pts2pxy(:,id_loc)*sig_loc;
% 
%     id_bin_y = mod(ibin,nbin)+1;
%     id_bin_x = ceil(ibin/nbin);
% 
%     id_grid_x = (id_bin_x*rbin-1) + (-((npxy-1)/2):((npxy-1)/2));
%     id_grid_y = (id_bin_y*rbin-1) + (-((npxy-1)/2):((npxy-1)/2));
%     id_grid_loc = (id_grid_x(:)-1)*ngrid + id_grid_y(:).';
% 
%     str2(id_grid_loc(:)) = str2(id_grid_loc(:)) + str_loc;
% end
% 
% toc;
% 
% pts2pxy = gs_kern(struct('r',r_loc),struct('r',pxy));
% 
% pts2grid = str2pxy\pts2pxy;
% str3 = zeros(ngrid*ngrid,1);
% 
% for ibin = 1:nbin*nbin
%     if idstart(ibin) == idstart(ibin+1), continue, end
% 
%     id_loc = idstart(ibin):(idstart(ibin+1)-1);
%     sig_loc = sigma_sort(id_loc);
% 
%     str_loc = pts2grid(:,id_loc)*sig_loc;
% 
%     id_bin_y = mod(ibin,nbin)+1;
%     id_bin_x = ceil(ibin/nbin);
% 
%     id_grid_x = (id_bin_x*rbin-1) + (-((npxy-1)/2):((npxy-1)/2));
%     id_grid_y = (id_bin_y*rbin-1) + (-((npxy-1)/2):((npxy-1)/2));
%     id_grid_loc = (id_grid_x(:)-1)*ngrid + id_grid_y(:).';
% 
%     str3(id_grid_loc(:)) = str3(id_grid_loc(:)) + str_loc;
% end
% 
% toc;
% 
% 
% 
% 
% 
% 
% figure(3);clf; t= tiledlayout('flow'); t.TileSpacing = 'tight';
% nexttile()
% my_scatter(rgrid,[],abs(str),'.')
% colorbar
% nexttile()
% my_scatter(rgrid,[],abs(str2),'.')
% colorbar
% nexttile()
% my_scatter(rgrid,[],log10(abs(str-str2)),'.')
% colorbar
% nexttile()
% my_scatter(rgrid,[],log10(abs(str3-str2)),'.')
% colorbar
% 
% 
% %% Matrix apply
% 
% kernvals = gs_kern(struct('r',[0;0;0]), struct('r',rgrid));
% 
% % fkern = kernel('lap','s');
% % kernvals = fkern.eval(struct('r',[0;0;0]), struct('r',rgrid));
% 
% kernvals = reshape(kernvals, ngrid,ngrid);
% 
% kernvalshift = fftshift(kernvals);
% % kernvalshift(1) = 0;
% 
% figure(4)
% subplot(1,2,1)
% my_scatter(rgrid,[],abs(kernvals(:)),'.')
% colorbar
% subplot(1,2,2)
% my_scatter(rgrid,[],abs(kernvalshift(:)),'.')
% colorbar
% 
% 
% tic;
% kern_hat = fft2(kernvalshift);
% str_hat = fft2(reshape(str,ngrid,ngrid));
% 
% u_hat = kern_hat.*str_hat;
% toc;
% 
% ugrid = ifft2(u_hat);
% 
% tic;
% utrue = zeros(ngrid,ngrid);
% str_sq = reshape(str,ngrid,ngrid);
% for i = 1:ngrid
%     for j = 1:ngrid
%         utrue(i,j) = sum(str_sq .* circshift(kernvalshift,[i-1,j-1]),"all");
%     end
% end
% toc;
% 
% ffterr = norm(utrue-ugrid)
% 
% figure(5);clf; t= tiledlayout('flow'); t.TileSpacing = 'tight';
% nexttile
% my_scatter(rgrid,[],real(ugrid(:))); colorbar
% 
% nexttile
% my_scatter(rgrid,[],real(utrue(:))); colorbar
% 
% % nexttile
% % my_scatter(rgrid,[],real(utrue(:)-u(:))); colorbar
% 
% nexttile
% my_scatter(rgrid,[],abs(str(:))); colorbar
% 
% 
% %% interp
% 
% u = zeros(S.npts,1);
% 
% 
% for ibin = 1:nbin*nbin
%     if idstart(ibin) == idstart(ibin+1), continue, end
% 
%     id_loc = idstart(ibin):(idstart(ibin+1)-1);
% 
%     id_bin_y = mod(ibin,nbin)+1;
%     id_bin_x = ceil(ibin/nbin);
% 
%     id_grid_x = (id_bin_x*rbin-1) + (-((npxy-1)/2):((npxy-1)/2));
%     id_grid_y = (id_bin_y*rbin-1) + (-((npxy-1)/2):((npxy-1)/2));
%     id_grid_loc = (id_grid_x(:)-1)*ngrid + id_grid_y(:).';
% 
%     u(id_loc) = u(id_loc) + pts2grid(:,id_loc).'*ugrid(id_grid_loc(:));
%     % str3(id_grid_loc(:)) = str3(id_grid_loc(:)) + str_loc;
% end
% 
% my_scatter(rsort,[],real(u))
% colorbar
% 

%% full apply setup

npxy = 9;
dx = 2*L/ngrid;
rpxy = 2 * (npxy-1)/2 * (dx);
tic;
[pxy,str2pxy,pxygrid] = setup_proxy(dx,npxy,rpxy,gs_kern,1e-8,zk_ext);
tproxy = toc

rbin = 3;
nbin = floor(ngrid/rbin);

tic;
[idstart,isort,rsort] = bin_pts(S.r,rbin,dx,nbin,L);


truebin = floor(1*nbin/6):ceil(5*nbin/6);
bin_ctrs = zeros(3,nbin*nbin);
for i = 1:nbin
    xcenter = rbin * (i-1) * dx - L;
    for j = 1:nbin
        ycenter =  rbin * (j-1) * dx - L;
        rcenter = [xcenter;ycenter;0];
        bin_ctrs(:,(i-1)*nbin + j)= rcenter;
    end
end
tbin = toc

tic;
r_loc = 0*rsort;
for ibin = 1:nbin*nbin
    if idstart(ibin) == idstart(ibin+1), continue, end

    id_loc = idstart(ibin):(idstart(ibin+1)-1);
    r_loc(:,id_loc) = rsort(:,id_loc) - bin_ctrs(:,ibin);
end

pts2pxy = gs_kern(struct('r',r_loc),struct('r',pxy));

pts2grid = str2pxy\pts2pxy;

Aspread = sparse(ngrid*ngrid,S.npts);
for id_bin_x = truebin
    for id_bin_y = truebin
    ibin = (id_bin_x-1)*nbin + id_bin_y;
    if idstart(ibin) == idstart(ibin+1), continue, end

    id_loc = idstart(ibin):(idstart(ibin+1)-1);

    id_grid_x = (id_bin_x-1)*rbin + (-((npxy-1)/2):((npxy-1)/2))+1;
    id_grid_y = (id_bin_y-1)*rbin + (-((npxy-1)/2):((npxy-1)/2))+1;

    id_grid_loc = (id_grid_x(:).'-1)*ngrid + id_grid_y(:);

    Aspread(id_grid_loc(:),id_loc) = Aspread(id_grid_loc(:),id_loc) + pts2grid(:,id_loc);
    end
end


tpre = toc


tic;
kernvals = gs_kern(struct('r',[0;0;0]), struct('r',rgrid));
kernvals = reshape(kernvals, ngrid,ngrid);

kernvalshift = fftshift(kernvals);
kern_hat = fft2(kernvalshift);
tkern = toc


%% full spread apply
t2 = tic;

tic;
str = Aspread*sigma(isort); str = full(str);
tspread = toc

tic;
str_hat = fft2(reshape(str,ngrid,ngrid));
u_hat = kern_hat.*str_hat;
ugrid = ifft2(u_hat);
tconv = toc

tic;
u = Aspread.'*ugrid(:);
tinterp = toc

tapply = toc(t2)
uplot = reshape(gs_kern(struct('r',S.r(:,1)), struct('r',rgrid))*sigma(1),ngrid,ngrid);
figure(3)
subplot(1,3,1)
my_scatter(rsort,[],u/norm(uplot,inf),'.'); colorbar

uplot = gs_kern(struct('r',S.r(:,1)),struct('r',rsort));
subplot(1,3,2)
my_scatter(rsort,[],uplot/norm(uplot,inf),'.');colorbar
subplot(1,3,3)
my_scatter(rsort,[],log10(abs(u-uplot)/norm(uplot,inf)),'.'); colorbar


%% redo setup with 2 points
% t1 = tic;

npxy = 9;
dx = 2*L/ngrid;
rpxy = 2 * (npxy-1)/2 * (dx);
% tic;
[pxy,str2pxy,pxygrid] = setup_proxy(dx,npxy,rpxy,gs_kern,1e-8,zk_ext);
% tproxy = toc

rbin = 4;
nbin = floor(ngrid/rbin);
S = []; 
% S.r = [[10;0;0],[0;-25;0]]; S.npts = size(S.r,2);
% S.r = [[12*dx;0;0],[0;-9*dx;0]]; S.npts = size(S.r,2); sigma = [1;-1];
% S.r = [[1.2*dx;5*dx;0]];S.npts = size(S.r,2); sigma = 1;
npts =20;% npts = 400;
S.r = 1.9*rad*[rand(2,npts)-0.5;zeros(1,npts)];S.npts = size(S.r,2); sigma = randn(npts,1);

% tic;
[idstart,isort,rsort] = bin_pts(S.r,rbin,dx,nbin,L);


truebin = floor(1*nbin/6):ceil(5*nbin/6);
bin_ctrs = zeros(3,nbin*nbin);
for i = 1:nbin
    xcenter = rbin * (i-1) * dx - L;
    for j = 1:nbin
        ycenter =  rbin * (j-1) * dx - L;
        rcenter = [xcenter;ycenter;0];
        bin_ctrs(:,(i-1)*nbin + j)= rcenter;
    end
end
% tbin = toc

% tic;
r_loc = 0*rsort;
for ibin = 1:nbin*nbin
    if idstart(ibin) == idstart(ibin+1), continue, end

    id_loc = idstart(ibin):(idstart(ibin+1)-1);
    r_loc(:,id_loc) = rsort(:,id_loc) - bin_ctrs(:,ibin);
end

pts2pxy = gs_kern(struct('r',r_loc),struct('r',pxy));

pts2grid = str2pxy\pts2pxy;

Aspread = sparse(ngrid*ngrid,S.npts);
for id_bin_x = truebin
    for id_bin_y = truebin
    ibin = (id_bin_x-1)*nbin + id_bin_y;
    if idstart(ibin) == idstart(ibin+1), continue, end

    id_loc = idstart(ibin):(idstart(ibin+1)-1);

    id_grid_x = (id_bin_x-1)*rbin + (-((npxy-1)/2):((npxy-1)/2))+1;
    id_grid_y = (id_bin_y-1)*rbin + (-((npxy-1)/2):((npxy-1)/2))+1;

    id_grid_loc = (id_grid_x(:).'-1)*ngrid + id_grid_y(:);

    Aspread(id_grid_loc(:),id_loc) = Aspread(id_grid_loc(:),id_loc) + pts2grid(:,id_loc);
    end
end

% tpre = toc


% tic;
kernvals = gs_kern(struct('r',[0;0;0]), struct('r',rgrid));
kernvals = reshape(kernvals, ngrid,ngrid);

kernvalshift = fftshift(kernvals);
kern_hat = fft2(kernvalshift);
% tkern = toc


%% full spread apply
t2 = tic;

% tic;
str = Aspread*sigma(isort); str = full(str);
% tspread = toc

% tic;
str_hat = fft2(reshape(str,ngrid,ngrid));
u_hat = kern_hat.*str_hat;
ugrid = ifft2(u_hat);
% tconv = toc

% tic;
u = Aspread.'*ugrid(:);
% tinterp = toc

% tapply = toc(t2)

% ttot = toc(t1)

% %%
% 
% u_slf = gs_kern(struct('r',[0;0;0]), struct('r',[0;0;0]));
% 
% u_fft = u - u_slf*sigma
% 
% tmp = gs_kern(struct('r',S.r), struct('r',S.r));
% u_dir = tmp(1,2)*sigma

%%

figure(6);clf; t= tiledlayout('flow'); t.TileSpacing = 'tight';
nexttile
% my_scatter(rgrid,[],real(ugrid(:)),'.'); colorbar
h = pcolor(xxgrid,yygrid,real(ugrid)); set(h,'edgecolor','none'); colorbar

% nexttile
% my_scatter(rsort,[],real(u(:))); colorbar

nexttile
% h = pcolor(xxgrid,yygrid,real(reshape(str,ngrid,ngrid))); set(h,'edgecolor','none'); colorbar
my_scatter(rgrid,[],real(str(:)),'o'); colorbar

hold on
my_scatter(S.r,'rx','linewidth',2)
id1 = round((S.r(1,:)+L)/(rbin*dx)+1);
id2 = round((S.r(2,:)+L)/(rbin*dx)+1);
idbin = (id1-1)*nbin + id2;
my_scatter(bin_ctrs(:,idbin),'ro','linewidth',2)

my_scatter(rgrid(:,id_grid_loc(:)),'r.','linewidth',2)
hold off
xlim([-6,6])
ylim([-6,6])

% nexttile
% my_scatter(rsort,[],real(sigma(isort))); colorbar

% nexttile
% my_scatter(S.r,[],real(sigma(:))); colorbar


uplot = reshape(gs_kern(struct('r',S.r), struct('r',rgrid))*sigma,ngrid,ngrid);
% uplot = reshape(gs_kern(struct('r',S.r([2,1,3])), struct('r',rgrid))*sigma,ngrid,ngrid);
nexttile
% my_scatter(rgrid,[],real(uplot(:)),'.'); colorbar
h = pcolor(xxgrid,yygrid,real(uplot)); set(h,'edgecolor','none'); colorbar

nexttile
% my_scatter(rgrid,[],real(uplot(:) - ugrid(:)),'.'); colorbar
h = pcolor(xxgrid,yygrid,log10(abs(uplot- ugrid))); set(h,'edgecolor','none'); colorbar

% %%
% % 
% 
% figure(5);clf; t= tiledlayout('flow'); t.TileSpacing = 'tight';
% nexttile
% h = pcolor(xxgrid,yygrid,real(ugrid)); set(h,'edgecolor','none'); colorbar
% 
% nexttile
% h = pcolor(xxgrid,yygrid,real(reshape(str,ngrid,ngrid))); set(h,'edgecolor','none'); colorbar
% 
% uplot = reshape(gs_kern(struct('r',S.r), struct('r',rgrid))*sigma,ngrid,ngrid);
% nexttile
% h = pcolor(xxgrid,yygrid,real(uplot)); set(h,'edgecolor','none'); colorbar
% 
% 
% src2pxy = gs_kern(struct('r',S.r), struct('r',pxy))*sigma;
% str_loc = str2pxy\src2pxy;
% 
% uplot2 = reshape(gs_kern(struct('r',pxygrid), struct('r',rgrid))*str_loc,ngrid,ngrid);
% nexttile
% h = pcolor(xxgrid,yygrid,real(uplot2)); set(h,'edgecolor','none'); colorbar
% 
% 
% nexttile
% h = pcolor(xxgrid,yygrid,log10(abs(ugrid- uplot2))); set(h,'edgecolor','none'); colorbar
% h = pcolor(xxgrid,yygrid,(real(ugrid- uplot2))); set(h,'edgecolor','none'); colorbar




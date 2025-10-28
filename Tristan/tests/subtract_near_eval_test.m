addpath(genpath('~/FastAlgs2/chunkie'))
addpath(genpath('~/FLAM'))
addpath(genpath('~/FastAlgs2/fmm3dbie_log_quad'))
% rng(123)

%% set parameters
alpha_int = 4;
beta_int = 1;
[rts_int,ejs_int] = surfwave.capillary.find_roots(alpha_int,beta_int,1);

alpha_ext = 1;
beta_ext = 100;
[rts_ext,ejs_ext] = surfwave.capillary.find_roots(alpha_ext,beta_ext,1);

zpars = [rts_ext;ejs_ext];

zk_int = rts_int(1);
zk_ext = rts_ext(1);

% gs_kern = @(s,t) surfwave.flex.kern(rts_ext,ejs_ext,s,t,'free plate eval second gs');
% gs_kern = @(s,t) surfwave.kern(rts_ext,ejs_ext,s,t,'gphi_s');
gs_kern_d = @(s,t) surfwave.kern(rts_ext,ejs_ext,s,t,'gs_d');
gs_kern_sp = @(s,t) surfwave.kern(rts_ext,ejs_ext,t,s,'gs_d').';
gs_kern_0 = @(s,t) surfwave.kern(rts_ext,ejs_ext,s,t,'gs_s');

%% Creating geometry and grid


rad = 1;
S = geometries.disk([rad,rad],[],[4 3 6],8);
% S = get_surfer(4,8,2*rad);

% ffun = @(t) starfish(t,5,0.,[],[],rad);
% cparams = []; cparams.maxchunklen = 4;
% chnkr = chunkerfunc(ffun,cparams);
% tic;
% S = get_surfer_chnkr(10,8,chnkr);
% tgeom = toc;

ngrid = 800;
L = 2.5*max(abs(S.r(:)));
xgrid = linspace(-L,L,ngrid+1); 
xgrid = xgrid(1:end-1);
% xgrid = xgrid(2:end);
[xxgrid, yygrid, zzgrid] = meshgrid(xgrid,xgrid,0);

rgrid = [xxgrid(:).'; yygrid(:).'; zzgrid(:).'];

ntarg = 4e4;
rtarg = [max(abs(S.r(:)));max(abs(S.r(:)));0] .* (rand(3,ntarg) - 0.5)*2;

srcinfo = []; srcinfo.r = S.r; srcinfo.n = [1;0;0] + 0*S.r;
targinfo = []; 
targinfo.r = rtarg; targinfo.r = S.r;
targinfo.n = [1;0;0] + 0*targinfo.r;

%% full apply setup

npxy = 9;
% npxy = 13;
dx = 2*L/ngrid;
rpxy = 2 * (npxy-1)/2 * (dx);
rbin = 2;


[Aspread_s,Aspread_t,Asubtract,kern_hat,rsort,isort,rsort_t,isort_t,times] = precom_eval_var(gs_kern_0,zk_ext,srcinfo,gs_kern_d,targinfo,gs_kern_0,gs_kern_d,L,dx,ngrid,rgrid,npxy,rpxy,rbin,1,1);
% [Aspread_s,Aspread_t,Asubtract,kern_hat,rsort,isort,rsort_t,isort_t,times] = precom_eval_var(gs_kern_0,zk_ext,srcinfo,gs_kern_0,targinfo,gs_kern_d,gs_kern_sp,L,dx,ngrid,rgrid,npxy,rpxy,rbin,0,1);

%% full spread apply
inz = randperm(S.npts,100);
% inz = 1;
sigma = zeros(S.npts,1); sigma(inz) = ones(length(inz),1);

u = apply_eval(sigma,ones(S.npts,1),Aspread_s,Aspread_t,Asubtract,kern_hat,isort,isort_t,ngrid, sparse(size(rtarg,2),S.npts));
u = u(isort_t);
% uplot = gs_kern(struct('r',S.r(:,inz)),struct('r',rsort_t))*sigma(inz);
uplot = gs_kern_d(struct('r',srcinfo.r(:,inz),'n',srcinfo.n(:,inz)),rsort_t)*sigma(inz);
% uplot = gs_kern_sp(struct('r',srcinfo.r(:,inz),'n',srcinfo.n(:,inz)),rsort_t)*sigma(inz);
% u = u*norm(uplot,inf)/norm(u,inf);

figure(5)
subplot(1,3,1)
my_scatter(rsort_t.r,[],real(u/norm(uplot,inf)),'.'); colorbar
hold on
my_scatter(S.r(:,inz),'x','linewidth',2)
hold off

% uplot = gs_kern(struct('r',S.r(:,1)),struct('r',rsort));
subplot(1,3,2)
my_scatter(rsort_t.r,[],real(uplot/norm(uplot,inf)),'.');colorbar
subplot(1,3,3)
my_scatter(rsort_t.r,[],log10(abs(u-uplot)/norm(uplot,inf)),'.'); colorbar
clim([-12,-7])
% hold on, plot(chnkr,'linewidth',2), hold off

tic;
str = Aspread_s*sigma(isort); str = full(str);
tspread = toc

tic;
str_hat = fft2(reshape(str,ngrid,ngrid));
u_hat = kern_hat.*str_hat;
ugrid = ifft2(u_hat);
tconv = toc

tic;
u2 = Aspread_t.'*ugrid(:);
tinterp = toc
figure(4)
subplot(1,2,1)
my_scatter(rsort_t.r,[],log10(abs(u-uplot)/norm(uplot,inf)),'.'); colorbar
% clim([-12,-7])

subplot(1,2,2)
my_scatter(rsort_t.r,[],log10(abs(u2-uplot)/norm(uplot,inf)),'.'); colorbar
% clim([-12,-7])




%%
sigma = exp(-vecnorm(S.r).^2/100).';
t1 = tic;
u = apply_eval(sigma,S.wts,Aspread_s,Aspread_t,Asubtract,kern_hat,isort,isort_t,ngrid, sparse(size(rtarg,2),S.npts));
toc(t1);


% function u = apply_mat(sigma,Aspread,Asubtract,kern_hat,isort,ngrid)
% tic;
% str = Aspread*sigma(isort); str = full(str);
% tspread = toc
% 
% tic;
% str_hat = fft2(reshape(str,ngrid,ngrid));
% u_hat = kern_hat.*str_hat;
% ugrid = ifft2(u_hat);
% tconv = toc
% 
% tic;
% u = Aspread.'*ugrid(:);
% tinterp = toc
% 
% tic;
% u = u + Asubtract*sigma(isort);
% tadd_sub = toc
% end

function u = apply_mat(sigma,Aspread,Asubtract,kern_hat,isort,ngrid)
    str = Aspread*sigma(isort); str = full(str);
    str_hat = fft2(reshape(str,ngrid,ngrid));
    u_hat = kern_hat.*str_hat;
    ugrid = ifft2(u_hat);
    u = Aspread.'*ugrid(:) + Asubtract*sigma(isort);
end

% function [Aspread,Asubtract,kern_hat,rsort,isort,times] = precom_apply(gkern,S,L,dx,ngrid,rgrid,npxy,rpxy,rbin,iftimes,tol)
% if nargin < nargin('precom_apply')-1
%     iftimes = 0;
% end
% 
% if nargin < nargin('precom_apply')
%     tol = 1e-10;
% end
% 
% tic;
% [pxy,pxy2str,pxygrid] = setup_proxy(dx,npxy,rpxy,gkern,tol);
% tproxy = toc;
% 
% 
% nbin = floor(ngrid/rbin);
% 
% tic;
% [idstart,isort,rsort] = bin_pts(S.r,rbin,dx,nbin,L);
% 
% 
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
% tbin = toc;
% 
% tic;
% r_loc = 0*rsort;
% for ibin = 1:nbin*nbin
%     if idstart(ibin) == idstart(ibin+1), continue, end
% 
%     id_loc = idstart(ibin):(idstart(ibin+1)-1);
%     r_loc(:,id_loc) = rsort(:,id_loc) - bin_ctrs(:,ibin);
% end
% 
% pts2pxy = gkern(struct('r',r_loc),struct('r',pxy));
% 
% pts2grid = pxy2str*pts2pxy;
% 
% Aspread = sparse(ngrid*ngrid,S.npts);
% for id_bin_x = truebin
%     for id_bin_y = truebin
%     ibin = (id_bin_x-1)*nbin + id_bin_y;
%     if idstart(ibin) == idstart(ibin+1), continue, end
% 
%     id_loc = idstart(ibin):(idstart(ibin+1)-1);
% 
%     id_grid_x = (id_bin_x-1)*rbin + (-((npxy-1)/2):((npxy-1)/2))+1;
%     id_grid_y = (id_bin_y-1)*rbin + (-((npxy-1)/2):((npxy-1)/2))+1;
% 
%     id_grid_loc = (id_grid_x(:).'-1)*ngrid + id_grid_y(:);
% 
%     Aspread(id_grid_loc(:),id_loc) = Aspread(id_grid_loc(:),id_loc) + pts2grid(:,id_loc);
%     end
% end
% 
% tpre = toc;
% 
% 
% tic;
% kernvals = gkern(struct('r',[0;0;0]), struct('r',rgrid));
% kernvals = reshape(kernvals, ngrid,ngrid);
% 
% kernvalshift = fftshift(kernvals);
% kern_hat = fft2(kernvalshift);
% tkern = toc;
% 
% nsub_rad = ceil(rpxy/dx);
% 
% nsub_bin = ceil((nsub_rad+(npxy-1)/2) / rbin);
% 
% nloc_grid = (nsub_bin*rbin + (npxy-1)/2);
% [rgrid_sub_x,rgrid_sub_y] = meshgrid((-nloc_grid:nloc_grid)*dx);
% rgrid_sub = [rgrid_sub_x(:).';rgrid_sub_y(:).';0*rgrid_sub_x(:).'];
% 
% [rbin_sub_x,rbin_sub_y] = meshgrid((-nsub_bin:nsub_bin)*rbin*dx);
% rbin_sub = [rbin_sub_x(:).';rbin_sub_y(:).';0*rbin_sub_x(:).'];
% 
% % figure(1);clf
% % my_scatter(rgrid_sub)
% % hold on
% % my_scatter(pxygrid,'.','linewidth',2)
% % my_scatter(pxy,'*','linewidth',2)
% % my_scatter(rbin_sub,'x','linewidth',2)
% % 
% % my_scatter(pxy+[nsub_bin*rbin*dx;0;0],'*','linewidth',2)
% % my_scatter(pxygrid+[nsub_bin*rbin*dx;0;0],'o','linewidth',2)
% % hold off
% 
% t1 = tic;
% loc_kern = gkern(struct('r',pxygrid),struct('r',rgrid_sub));
% % Asubtract = sparse(S.npts, S.npts);
% % for id_bin_x = truebin
% %     for id_bin_y = truebin
% %         ibin = (id_bin_x-1)*nbin + id_bin_y;
% %         if idstart(ibin) == idstart(ibin+1), continue, end
% % 
% %         id_loc = idstart(ibin):(idstart(ibin+1)-1);
% %         sig2loc = loc_kern*pts2grid(:,id_loc);
% % 
% %         rsrc = rsort(:,id_loc);
% % 
% %         ibins = (id_bin_x+(-nsub_bin:nsub_bin)-1)*nbin + id_bin_y+(-nsub_bin:nsub_bin).';
% %         ibins = ibins(:);
% %         idstart_loc = zeros(1,length(ibins)+1);
% %         idstart_loc(1) = 1;
% %         id_targs = [];
% %         for i = 1:length(ibins)
% %             idstart_loc(i+1) = idstart_loc(i) + idstart(ibins(i)+1)-idstart(ibins(i));
% %             id_targs = [id_targs,idstart(ibins(i)):(idstart(ibins(i)+1)-1)];
% %         end
% % 
% %         sig2u = gkern(struct('r',rsrc),struct('r',rsort(:,id_targs)));
% % 
% %         for i = -nsub_bin:nsub_bin
% %             for j = -nsub_bin:nsub_bin
% % 
% %             ibin = (id_bin_x+i-1)*nbin + id_bin_y+j;
% %             id_loc_targ = idstart(ibin):(idstart(ibin+1)-1);
% % 
% %             id_grid_x = i*rbin + (-((npxy-1)/2):((npxy-1)/2))+1+nloc_grid;
% %             id_grid_y = j*rbin + (-((npxy-1)/2):((npxy-1)/2))+1+nloc_grid;
% % 
% %             id_grid_loc = (id_grid_x(:).'-1)*(2*nloc_grid+1) + id_grid_y(:);
% % 
% %             sig2targ = pts2grid(:,id_loc_targ).' * sig2loc(id_grid_loc(:),:);
% %             % sig2uloc = gs_kern(struct('r',rsrc),struct('r',rsort(:,id_loc_targ)));
% % 
% %             ibinloc = (i+nsub_bin+1-1)*(2*nsub_bin+1) + j+nsub_bin+1;
% %             sig2uloc = sig2u(idstart_loc(ibinloc):(idstart_loc(ibinloc+1)-1),:);
% %             Asubtract(id_loc_targ,id_loc) = sig2uloc - sig2targ;
% % 
% %             end
% %         end
% %     end
% % end
% % toc;
% % Asubtract2 = Asubtract;
% Asubtract = sparse(S.npts, S.npts);
% tic;
% id_startbin = ones(1,nbin*nbin); id_startbin(1) = 1;
% rnear = [];
% for id_bin_x = truebin
%     for id_bin_y = truebin
%         ibin = (id_bin_x-1)*nbin + id_bin_y;
%         if idstart(ibin) == idstart(ibin+1),id_startbin(ibin+1:end) = id_startbin(ibin);continue, end
% 
%         id_loc = idstart(ibin):(idstart(ibin+1)-1);
%         % sig2loc = loc_kern*pts2grid(:,id_loc);
% 
%         rsrc = rsort(:,id_loc);
% 
%         ibins = (id_bin_x+(-nsub_bin:nsub_bin)-1)*nbin + id_bin_y+(-nsub_bin:nsub_bin).';
%         ibins = ibins(:);
%         idstart_loc = zeros(1,length(ibins)+1);
%         idstart_loc(1) = 1;
%         id_targs = [];
%         for i = 1:length(ibins)
%             idstart_loc(i+1) = idstart_loc(i) + idstart(ibins(i)+1)-idstart(ibins(i));
%             id_targs = [id_targs,idstart(ibins(i)):(idstart(ibins(i)+1)-1)];
%         end
%         rnear_loc = rsort(:,id_targs) - reshape(rsrc,3,1,[]);
%         % rnear_loc = reshape(rsort(:,id_targs),3,1,[]) - rsrc;
%         rnear = [rnear, rnear_loc(:,:)];
%         id_startbin(ibin+1) = id_startbin(ibin) + size(rnear_loc(:,:),2);
% 
% 
%     end
% end
% toc;
% tic;
%     sig2u = gkern(struct('r',[0;0;0]),struct('r',rnear));
% toc;
% tic;
% pts2subpts = loc_kern*pts2grid;
% for id_bin_x = truebin
%     for id_bin_y = truebin
%         ibin = (id_bin_x-1)*nbin + id_bin_y;
%         if idstart(ibin) == idstart(ibin+1), continue, end
% 
%         id_loc = idstart(ibin):(idstart(ibin+1)-1);
%         % sig2loc = loc_kern*pts2grid(:,id_loc);
%         sig2loc = pts2subpts(:,id_loc);
%         % a = pts2grid(id_loc,:).';
% 
%         % rsrc = rsort(:,id_loc);
% 
%         ibins = (id_bin_x+(-nsub_bin:nsub_bin)-1)*nbin + id_bin_y+(-nsub_bin:nsub_bin).';
%         ibins = ibins(:);
%         idstart_loc = zeros(1,length(ibins)+1);
%         idstart_loc(1) = 1;
%         for i = 1:length(ibins)
%             idstart_loc(i+1) = idstart_loc(i) + idstart(ibins(i)+1)-idstart(ibins(i));
% 
%         end
% 
%         sig2ubin = sig2u(id_startbin(ibin):(id_startbin(ibin+1)-1));
%         sig2ubin = reshape(sig2ubin,[],length(id_loc));
%         for i = -nsub_bin:nsub_bin
%             for j = -nsub_bin:nsub_bin
% 
%             ibin = (id_bin_x+i-1)*nbin + id_bin_y+j;
%             id_loc_targ = idstart(ibin):(idstart(ibin+1)-1);
% 
%             id_grid_x = i*rbin + (-((npxy-1)/2):((npxy-1)/2))+1+nloc_grid;
%             id_grid_y = j*rbin + (-((npxy-1)/2):((npxy-1)/2))+1+nloc_grid;
% 
%             id_grid_loc = (id_grid_x(:).'-1)*(2*nloc_grid+1) + id_grid_y(:);
% 
%             sig2targ = pts2grid(:,id_loc_targ).' * sig2loc(id_grid_loc(:),:);
%             % sig2uloc = gs_kern(struct('r',rsrc),struct('r',rsort(:,id_loc_targ)));
% 
%             ibinloc = (i+nsub_bin+1-1)*(2*nsub_bin+1) + j+nsub_bin+1;
%             sig2uloc = sig2ubin(idstart_loc(ibinloc):(idstart_loc(ibinloc+1)-1),:);
%             Asubtract(id_loc_targ,id_loc) = sig2uloc - sig2targ;
%             end
%         end
%     end
% end
% toc;
% tsub_pre = toc(t1);
% 
% 
% times = [tsub_pre,tpre,tproxy,tbin,tkern];
% 
% if iftimes
%     % %%
%     % fprintf('add and subtract op =  %.2e s\n',tsub_pre)
%     % fprintf('spread op = %.2e s\n',tpre)
%     % fprintf('proxy = %.2e s\n',tproxy)
%     % fprintf('binning = %.2e s\n',tbin)
%     % fprintf('fft kern =  %.2e s\n',tkern)
% 
%     tproxy
%     tbin
%     tkern
%     tpre
%     tsub_pre
%     % %%
% end
% 
% end


function S = get_surfer(hmax,norder,rad)
% gm = fegeometry(@lshapeg);
t = pi/12:pi/12:2*pi;
pgon = polyshape({rad*[-0.5 -0.5 0.5 0.5], rad*0.25*cos(t)}, ...
                 {rad*[0.5 -0.5 -0.5 0.5], rad*0.25*sin(t)});
tr = triangulation(pgon);
gm = fegeometry(tr);
% gm = fegeometry(@lshapeg);
gm = generateMesh(gm,GeometricOrder='linear',Hmax=hmax);
% 
npatches = size(gm.Mesh.Elements,2);

triaskel = zeros(3,3,npatches);
for i = 1:npatches
    triaskel(1:2,:,i) = gm.Mesh.Nodes(:,gm.Mesh.Elements(:,i));
end

% [npols, uvs] =  vioreanu_simplex_quad(norder);
[uvs]=koorn.rv_nodes(norder); npols = size(uvs,2);
[srcvals,npts] = triangles_to_points(npatches, triaskel,npols, uvs);

iptype = ones(npatches,1);
S = surfer(npatches, norder, srcvals, iptype);

% pdemesh(gm), hold on, my_scatter(S.r,'.'), hold off
end




function [srcvals,npts] = triangles_to_points(npatches, triaskel,npols, uvs)
npts = npatches*npols;
srcvals = zeros(12,npts);

for ipatch = 1:npatches
  for i = 1:npols
    u = uvs(1,i);
    v = uvs(2,i);
    
    ipt = (ipatch-1)*npols + i;
    
    [srcvals(1:3,ipt),srcvals(4:9,ipt)] =  xtri_flat_eval(ipatch,u,v,triaskel);

    srcvals(10:12,ipt) =  cross(srcvals(4:6,ipt),srcvals(7:9,ipt));

    ds = sqrt(srcvals(10,ipt)^2 + srcvals(11,ipt)^2 + srcvals(12,ipt)^2);
    srcvals(10,ipt) = srcvals(10,ipt)/ds;
    srcvals(11,ipt) = srcvals(11,ipt)/ds;
    srcvals(12,ipt) = srcvals(12,ipt)/ds;

  end
  ipts = (ipatch-1)*npols + (1:npols);

end

end





function [xyz,dxyzduv] = xtri_flat_eval(itri, u, v, triainfo)

  xyz = zeros(3,1);
  dxyzduv = zeros(3,2);

  %
  % project the triangle itri in triainfo onto the sphere
  %
  %    Input:
  % itri - triangle number to map
  % u,v - local uv coordinates on triangle itri
  % triainfo - flat skeleton triangle info
  %
  %    Output:
  % xyz - point on the sphere
  % dxyzduv - first derivative information
  %
  %s

  x0=triainfo(1:3,1,itri);

  x1=triainfo(1:3,2,itri);

  x2=triainfo(1:3,3,itri);

  %
  % ... process the geometry, return the point location on the sphere
  % and the derivatives with respect to u and v
  %
  xyz=x0+u*(x1-x0)+v*(x2-x0);

  dxyzduv(:,1) = x1-x0;

  dxyzduv(:,2) = x2-x0; 
  
  dxyzduv = dxyzduv(:);
end

function S = get_surfer_chnkr(hmax,norder,chnkr)

[rend,tauend] = chunkends(chnkr);
verts = reshape(rend(1:2,1,:),2,[]);

pgon = polyshape(verts(1,:), verts(2,:));
tr = triangulation(pgon);
gm = fegeometry(tr);
gm = addVertex(gm,Coordinates=tr.Points);
% gm = fegeometry(@lshapeg);
gm = generateMesh(gm,GeometricOrder='linear',Hmax=hmax);
% 
npatches = size(gm.Mesh.Elements,2);

triaskel = zeros(3,3,npatches);
for i = 1:npatches
    triaskel(1:2,:,i) = gm.Mesh.Nodes(:,gm.Mesh.Elements(:,i));
end

% [npols, uvs] =  vioreanu_simplex_quad(norder);
[uvs]=koorn.rv_nodes(norder); npols = size(uvs,2);
[srcvals,npts] = triangles_to_points_chnkr(npatches, triaskel,npols, uvs,chnkr);
% [srcvals,npts] = triangles_to_points(npatches, triaskel,npols, uvs);

iptype = ones(npatches,1);
S = surfer(npatches, norder, srcvals, iptype);

pdemesh(gm), hold on, my_scatter(S.r,'.'),plot(chnkr), hold off
end

function [srcvals,npts] = triangles_to_points_chnkr(npatches, triaskel,npols, uvs,chnkr)
npts = npatches*npols;
srcvals = zeros(12,npts);

for ipatch = 1:npatches
  for i = 1:npols
    u = uvs(1,i);
    v = uvs(2,i);
    
    ipt = (ipatch-1)*npols + i;
    
    [srcvals(1:3,ipt),srcvals(4:9,ipt)] =  xtri_flat_eval(ipatch,u,v,triaskel);

    srcvals(10:12,ipt) =  cross(srcvals(4:6,ipt),srcvals(7:9,ipt));

    ds = sqrt(srcvals(10,ipt)^2 + srcvals(11,ipt)^2 + srcvals(12,ipt)^2);
    srcvals(10,ipt) = srcvals(10,ipt)/ds;
    srcvals(11,ipt) = srcvals(11,ipt)/ds;
    srcvals(12,ipt) = srcvals(12,ipt)/ds;

  end

end

[rend,~] = chunkends(chnkr);

for i = 1:chnkr.nch

    id1s = any(vecnorm(triaskel(1:2,:,:) - rend(1:2,1,i)) < 1e-12);
    id2s = any(vecnorm(triaskel(1:2,:,:) - rend(1:2,2,i)) < 1e-12);


    ipatch = find(all([id1s(:).'; id2s(:).']));
    
    triloc = triaskel(:,:,ipatch);

    [~,id1] = min(vecnorm(triloc(1:2,:) - rend(1:2,1,i)));
    [~,id2] = min(vecnorm(triloc(1:2,:) - rend(1:2,2,i)));

    id3 = setdiff(1:3,[id1,id2]);

    triloc = triloc(:, [id1,id3,id2]);

    rloc = chnkr.r(:,:,i);
    dloc = chnkr.d(:,:,i);
    if size(rloc,1) == 2
        rloc = [rloc;zeros(1,16)];
        dloc = [dloc;zeros(1,16)];
    end

    for j = 1:npols
        u = uvs(1,j);
        v = uvs(2,j);
        
        ipt = (ipatch-1)*npols + j;
        
        [srcvals(1:3,ipt),srcvals(4:9,ipt)] =  xtri_param_eval(triloc,u,v,rloc,dloc);
    
        srcvals(10:12,ipt) =  cross(srcvals(4:6,ipt),srcvals(7:9,ipt));
    
        ds = sqrt(srcvals(10,ipt)^2 + srcvals(11,ipt)^2 + srcvals(12,ipt)^2);
        srcvals(10,ipt) = srcvals(10,ipt)/ds;
        srcvals(11,ipt) = srcvals(11,ipt)/ds;
        srcvals(12,ipt) = srcvals(12,ipt)/ds;
    
    end
end



end


function [xyz,dxyzduv] = xtri_param_eval(triainfo, u, v,rs,ds)

  xyz = zeros(3,1);
  dxyzduv = zeros(3,2);


  %
  % project the triangle itri in triainfo onto the sphere
  %
  %    Input:
  % itri - triangle number to map
  % u,v - local uv coordinates on triangle itri
  % triainfo - flat skeleton triangle info
  %
  %    Output:
  % xyz - point on the sphere
  % dxyzduv - first derivative information
  %
  %s

  x0=triainfo(1:3,1);

  x1=triainfo(1:3,2);

  x2=triainfo(1:3,3);

    % rs = rs - x1;
    k = size(rs,2);
    [~,~,val2coef] = lege.exps(k);
    
    cri = val2coef*(rs.');
    cdi = val2coef*(ds.');

    s = (1-u);
    t = 1 - v/(1-u);
    
    legs = lege.pols(2*t-1,k-1); legs = legs.';
    
    rv = legs*cri; rv = rv.';
    dv = legs*cdi; dv = -dv.';



  %
  % ... process the geometry, return the point location on the sphere
  % and the derivatives with respect to u and v
  %
  % xyz=x0+u*(x1-x0)+rv;
  % 
  % dxyzduv(:,1) = x1-x0;

    xyz=(rv-x1)*s + x1;

  % dxyzduv(:,1) = (-(rv-x1) + t*dv);
  dxyzduv(:,1) = (-(rv-x1) + v/(1-u)*dv);

  dxyzduv(:,2) = -dv; 
  
  dxyzduv = dxyzduv(:);

end

% fmm3dbie integration demo
% Requires fmm3dbie: https://github.com/fastalgorithms/fmm3dbie
%
% Solve Dirichlet scattering problems with several inclusions
% Takes about 10 minutes

% Need to load planewave function
addpath("utils")

% ambient wavenumber
zk = 10;
kvec = zk*[1;0;0];

%% Make a field of scatters

ntry = 1000;
rad = 1;
% make interior boundaries with random locations
surfs = [];
L = 4;
ctrs = zeros(3,1);
nscat = 5;
for i  = 1:nscat
    % each boundary is sphere
    surfi_i = geometries.sphere(rad,3,ctrs(:,i),6);
    surfs = [surfs,surfi_i];

    % try to find location for next chunker
    for j = 1:ntry
        tmp = L*(2*rand(3,1)-1);
        rmin = min(vecnorm(tmp - ctrs));
        if (rmin > rad * 3); break; end
    end
    if j == ntry; error('Could not place next boundary'); end
    ctrs = [ctrs,tmp];
end
ctrs = ctrs(:,1:end-1);
S = merge(surfs);
npt_int = S.npts;

fprintf('Geometry generated\n')

%% Make system
% define system kernel

skern = @(s,t) helm3d.kern(zk, s, t, 's');
dkern = @(s,t) helm3d.kern(zk, s, t, 'd');

% whether to run the accuracy test
ifacc = 0;

% define boundary data
if ifacc
    src = struct('r',[0;0;0]);
    rhs = -skern(src,S);
else
    rhs = -planewave(kvec,S.r(:,:));
    rhs = rhs(:);
end

figure(3); clf;
plot(S,real(-rhs));colorbar
axis equal
drawnow()
%% Compute example field

L = 1.1*max(vecnorm(S.r(:,:)));
x1 = linspace(-L,L,100);
[xx,yy] = meshgrid(x1,x1);
targs = [xx(:).'; yy(:).';0*xx(:).'];
ntargs = size(targs,2);

% identify points in computational domain
in = any(((targs(1,:)-ctrs(1,:).').^2+(targs(2,:)-ctrs(2,:).').^2+(targs(3,:)-ctrs(3,:).').^2)<rad);
out = ~in;

targout = [];
targout.r = targs(:,out);

% get incoming solution
uin = zeros(size(xx));
if ifacc
    uin(out) = skern(src,targout);
else
    uin(out) = planewave(kvec(:),targout.r);
end


%% Precompuation

eps = 1e-6;

t1 = tic;
Q = helm3d.dirichlet.get_quadrature_correction(S,eps,zk,[0,1],S);
novers = get_oversampling_parameters(S,Q,eps);
[S_over,interp_mat] = oversample(S,novers);
cors = d_cor_mat(S, S_over, interp_mat, zk, Q);

% determine a good spreading box size, 
[grid_info, ~] = get_grid(skern, S, S, eps);
proxy_opts = [];
proxy_opts.halfside = grid_info.dx*grid_info.nbinpts;
% setup real grid
[grid_info, proxy_info] = get_grid(skern, S_over, targout, eps, [], proxy_opts);
% get spreading operators
[A_spread_s, sort_info_S_over]= get_spread(skern, dkern, S_over, ...
    grid_info, proxy_info, {'r','n'});
[A_spread_surf, sort_info_S]= get_spread(skern, skern, S, ...
    grid_info, proxy_info);
[A_spread_t, sort_info_t] = get_spread(skern, skern, targout, ...
    grid_info, proxy_info);
% build corrections
A_addsub_self = get_addsub(skern, dkern, S_over, S, ...
    grid_info, proxy_info, sort_info_S_over, sort_info_S, A_spread_s, A_spread_surf);
A_addsub_eval = get_addsub(skern, dkern, S_over, targout, ...
    grid_info, proxy_info, sort_info_S_over, sort_info_t, A_spread_s, A_spread_t);
% get DFT of kernel
skern_hat = get_kernhat(skern,grid_info);

% %% compare with dense
% tic;
% Q = helm3d.dirichlet.get_quadrature_correction(S,eps,zk,[0,1],S,struct('format','sparse'));
% % %%
% P = zeros(S.npts,1);
% Amat = helm3d.dirichlet.matgen(1:S.npts,1:S.npts,S,[zk,0,1],P,Q);
% tdens = toc
% submat = helm3d.kern(zk, S, S, 'd').*S.wts.';
% submat(isnan(submat)) = 0;
% Bmat = cors + submat + 0.5*eye(size(cors));
% Bmat = submat;
% idx = abs(Q.spmat) ~=0;
% Bmat(idx) = Q.spmat(idx); 
% Bmat = Bmat + 0.5*eye(size(cors));

cors = cors + (A_addsub_self.*S_over.wts(:).')*interp_mat + 0.5*speye(size(cors));
A_spread_s = (A_spread_s.*S_over.wts(:).') * interp_mat;
A_addsub_eval = (A_addsub_eval.*S_over.wts(:).') * interp_mat;
tprecom = toc(t1)

%%
sys_app = @(dens) pcfft_apply(dens,A_spread_s,A_spread_surf,cors,skern_hat);

tic;
% solve
sol = gmres(sys_app,rhs,[],eps,1000);
tsolve = toc

%% compute utot on slice of the 2D domain for plotting

% get solution
tic;
uscat = (NaN+NaN*1i)*zeros(size(xx));
Q = helm3d.dirichlet.get_quadrature_correction(S,eps,zk,[0,1],targout);
corsplot = d_cor_mat(S, S_over, interp_mat, zk, Q);
uscat(out) = pcfft_apply(sol, A_spread_s, A_spread_t,corsplot+A_addsub_eval,skern_hat);
tplot = toc

utot = uin + uscat;

%% make plots
umax = max(abs(utot(:))); 
figure(2);clf
h = pcolor(xx,yy,imag(utot)); set(h,'EdgeColor','none'); colorbar
colormap(redblue); %clim(0.4*[-umax,umax]);
hold on 
h = plot(S);
h.FaceColor = '#999';
axis equal


function cormat = d_cor_mat(S, S_over, interp_mat, zk, Q)

Q2 = Q;
Q2.wnear = Q.wnear;

    ixyzs = S.ixyzs(:);
    ixyzs_over = S_over.ixyzs(:);
    npatches = S.npatches;

    istarts = Q.row_ptr(1:end-1);
    iends = Q.row_ptr(2:end)-1;
    isrcinds = cell(npatches,1);
    isrcinds_over = cell(npatches,1);
    for i=1:npatches
        isrcinds{i} = ixyzs(i):ixyzs(i+1)-1;
        isrcinds_over{i} = ixyzs_over(i):ixyzs_over(i+1)-1;
    end

    for i=1:size(Q.targinfo.r(:,:),2)
        iinds = horzcat(isrcinds{Q.col_ind(istarts(i):iends(i))});
        iinds_over = horzcat(isrcinds_over{Q.col_ind(istarts(i):iends(i))});
        srcinfo = [];
        srcinfo.r = S_over.r(:,iinds_over);
        srcinfo.n = S_over.n(:,iinds_over);
        targinfo = [];
        targinfo.r = Q.targinfo.r(:,i);
        % targinfo.n = Q.n(:,i);
        submat = helm3d.kern(zk, srcinfo, targinfo, 'd');
        submat(isnan(submat)) = 0;
        Q2.wnear(Q.iquad(istarts(i)):Q.iquad(iends(i)+1)-1) = Q2.wnear(Q.iquad(istarts(i)):Q.iquad(iends(i)+1)-1) - interp_mat(iinds_over,iinds).' * (submat.'.*S_over.wts(iinds_over));
    end

cormat = conv_rsc_to_spmat(S, Q2.row_ptr, Q2.col_ind, Q2.wnear);

end
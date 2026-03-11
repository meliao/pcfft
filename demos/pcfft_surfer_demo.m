% Chunkie integration demo
% Requires chunkIE: https://github.com/fastalgorithms/chunkie
%
% Solve Dirichlet scattering problems with many inclusions

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
ctrs = L*rand(3,1);
nscat = 2;
for i  = 1:nscat
    % each boundary is sphere
    surfi_i = geometries.sphere(rad,4,ctrs(:,i));
    surfs = [surfs,surfi_i];

    % try to find location for next chunker
    for j = 1:ntry
        tmp = L*rand(3,1);
        rmin = min(vecnorm(tmp - ctrs));
        if (rmin > rad * 3); break; end
    end
    if j == ntry; error('Could not place next boundary'); end
    ctrs = [ctrs,tmp];
end
S = merge(surfs);
npt_int = S.npts;

fprintf('Geometry generated\n')
figure(3); clf;
plot(S)
axis equal

%% Make system
% define system kernel

skern = @(s,t) helm3d.kern(zk, s, t, 's');
dkern = @(s,t) helm3d.kern(zk, s, t, 'd');
% spkern = kernel('helm', 'sp', zk);

% define boundary data
rhs = -planewave(kvec,S.r(:,:));
rhs = rhs(:);
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
uin(out) = planewave(kvec(:),targs(:,out));

%% Precompuation

eps = 1e-6;


t1 = tic;

% setup grid
[grid_info, proxy_info] = get_grid(skern, S, targout, eps);
% get spreading operators
[A_spread_s, sort_info_c]= get_spread(skern, dkern, S, ...
    grid_info, proxy_info, {'r','n'});
[A_spread_surf, ~]= get_spread(skern, skern, S, ...
    grid_info, proxy_info);
[A_spread_t, sort_info_t] = get_spread(skern, skern, targout, ...
    grid_info, proxy_info);
% build corrections
A_addsub_c = get_addsub(skern, dkern, S, S, ...
    grid_info, proxy_info, sort_info_c, sort_info_c, A_spread_s, A_spread_surf);
A_addsub_eval = get_addsub(skern, dkern, S, targout, ...
    grid_info, proxy_info, sort_info_c, sort_info_t, A_spread_s, A_spread_t);
% get DFT of kernel
skern_hat = get_kernhat(skern,grid_info);


opts = [];
opts.format = 'sparse';
% not actually corrections
cors = helm3d.dirichlet.get_quadrature_correction(S,eps,zk,[0,1],S,opts);
cors = cors.spmat + A_addsub_c.*S.wts(:).' + 0.5*speye(size(cors.spmat));
A_spread_s = A_spread_s.*S.wts(:).';
A_addsub_eval = A_addsub_eval.*S.wts(:).';
tprecom = toc(t1)

sys_app = @(dens) pcfft_apply(dens,A_spread_s,A_spread_surf,cors,skern_hat);
% return

tic;
% solve
sol = gmres(sys_app,rhs,[],eps,1000);
tsolve2 = toc

%% compute utot on the 2D domain for plotting

% get solution
tic;
uscat = (NaN+NaN*1i)*zeros(size(xx));
cors = helm3d.dirichlet.get_quadrature_correction(S,eps,zk,targout,[0,1],opts);
uscat(out) = pcfft_apply(sol, A_spread_s, A_spread_t,cors.spmat+A_addsub_eval,skern_hat);
tplot = toc

utot = uin + uscat;

%% make plots
umax = max(abs(utot(:))); 
figure(2);clf
h = pcolor(xx,yy,imag(utot)); set(h,'EdgeColor','none'); colorbar
colormap(redblue); clim(0.4*[-umax,umax]);
hold on 
plot(surfs,'k')
axis equal
% title('$u^{\textrm{tot}}$','Interpreter','latex','FontSize',12)

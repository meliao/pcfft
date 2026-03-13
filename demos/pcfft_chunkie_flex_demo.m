%Chunkie integration demo2
% Requires chunkIE: https://github.com/fastalgorithms/chunkie
%
% Solve flexural scattering problems with many inclusions

addpath("utils")

% planewave direction
phi = 0;
% ambient wavenumber
zk = 10;
kvec = zk*[cos(phi);sin(phi)];

%% Make a field of scatters

% geometry preferences
cparams = [];
% cparams.maxchunklen = min(4.0/max(zk),0.125);
cparams.maxchunklen = min(4.0/max(zk));
pref = []; 
pref.k = 16;
narms =5;
amp = 0.25;
scale = .6;
rad = scale*(amp+1);
ntry = 1000;

% make interior boundaries with random locations
chnkr = [];
L = 8;
% L = 2;
theta = 2*pi*rand();
ctrs = L*rand()*[cos(theta);sin(theta)];
n_pts = [];
nscat = 15;
for i  = 1:nscat
    % each boundary is rotated starfish with a random number of arms
    phi = 2*pi*rand();
    narms = randi([3,6]);
    chnkr_i = chunkerfunc(@(t) starfish(t,narms,amp,ctrs(:,i),phi,scale), ...
        cparams,pref); 
    % track number of points
    n_pts(i) = chnkr_i.npt;
    chnkr = [chnkr,chnkr_i];

    % try to find location for next chunker
    for j = 1:ntry
        theta = 2*pi*rand();
        tmp = L*rand()*[cos(theta);sin(theta)];
        rmin = min(vecnorm(tmp - ctrs));
        if (rmin > rad * 3); break; end
    end
    if j == ntry; error('Could not place next boundary'); end
    ctrs = [ctrs,tmp];
end
ctrs = ctrs(:,1:end-1);
chnkr = merge(chnkr);
npt_int = chnkr.npt;

fprintf('Geometry generated\n')
figure(3); clf;
plot(chnkr)
axis equal

%% Make system
% define system kernel

skern = @(s,t) chnk.flex2d.kern(zk, s, t, 's');   
ekern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');   
bkern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_bcs');   
fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate');

% define boundary data
[r1, grad] = planewave(kvec, chnkr.r);

nx = chnkr.n(1,:); 
ny = chnkr.n(2,:);

normalderiv = grad(:, 1).*(nx.')+ grad(:, 2).*(ny.');                                % Dirichlet and Neumann BC(Clamped BC)                         

firstbc = -r1;
secondbc = -normalderiv;

rhs = zeros(2*chnkr.npt, 1); rhs(1:2:end) = firstbc ; rhs(2:2:end) = secondbc;

% src = []; src.r = [-10;0];
% rhs = -bkern(src,chnkr);

%% Compute example field

L = 1.1*max(vecnorm(chnkr.r(:,:)));
x1 = linspace(-L,L,300);
[xx,yy] = meshgrid(x1,x1);
targs = [xx(:).'; yy(:).'];
ntargs = size(targs,2);

% identify points in computational domain
in = chunkerinterior(chnkr,{x1,x1});
out = ~in;

targout = [];
targout.r = targs(:,out);
% get incoming solution
uin = zeros(size(xx));
uin(out) = planewave(kvec(:),targs(:,out));
% uin(out) =skern(src,targout);

%% Precompuation

eps = 1e-6;

srcs = [];
srcs.r = chnkr.r(:,:);
srcs.n = chnkr.n(:,:);

t1 = tic;

proxy_opts = [];
% as this is a fourth order problem, we proxy against values and first
% derivatives
proxy_opts.proxy_der = 1;

% setup grid
[grid_info, proxy_info] = get_grid(skern, chnkr, targout, eps,[],proxy_opts);
% get spreading operators
[A_spread_s, sort_info_c]= get_spread(skern, ekern, chnkr, ...
    grid_info, proxy_info, {'r','n','d'});
[A_spread_c, ~]= get_spread(skern, @(s,t) bkern(t,s).', chnkr, ...
    grid_info, proxy_info, {'r','n','d'});
[A_spread_t, sort_info_t] = get_spread(skern, skern, targout, ...
    grid_info, proxy_info);
% build corrections
A_addsub_c = get_addsub(skern, fkern, chnkr, chnkr, ...
    grid_info, proxy_info, sort_info_c, sort_info_c, A_spread_s, A_spread_c);
A_addsub_eval = get_addsub(skern, ekern, chnkr, targout, ...
    grid_info, proxy_info, sort_info_c, sort_info_t, A_spread_s, A_spread_t);
% get DFT of kernel
skern_hat = get_kernhat(skern,grid_info);


kappa = signed_curvature(chnkr);
kappa = kappa(:);
wts = repmat(chnkr.wts(:).',2,1);

opts = [];
opts.sing = 'log';
opts.corrections = 1;

% a = fkern(chnkr,chnkr).*wts(:).';
% a(1:2*size(a,1)+2:end) = 0;
% a(2:2*size(a,1)+2:end) = 0;
% a(1+size(a,1):2*size(a,1)+2:end) = 0;
% a(2+size(a,1):2*size(a,1)+2:end) = 0;

start = tic;
cors = chunkermat(chnkr,fkern, opts);
cors = cors - 0.5*speye(2*chnkr.npt) + A_addsub_c.*wts(:).';
% cors = cors - 0.5*speye(2*chnkr.npt) + a;
cors(2:2:end,1:2:end) = cors(2:2:end,1:2:end) + kappa.*speye(chnkr.npt);

A_spread_s = A_spread_s.*wts(:).';
A_addsub_eval = A_addsub_eval.*wts(:).';
tprecom = toc(t1)

sys_app = @(dens) pcfft_apply(dens,A_spread_s,A_spread_c,cors,skern_hat);

tic;
% solve
sol = gmres(sys_app,rhs,[],eps,1000);
tsolve2 = toc

% opts.corrections = 0;
% start = tic;
% sys = chunkermat(chnkr,fkern, opts);
% sys = sys - 0.5*eye(2*chnkr.npt);
% sys(2:2:end,1:2:end) = sys(2:2:end,1:2:end) + kappa.*eye(chnkr.npt);
% sol2 = sys\rhs;

%%

% get solution
tic;
uscat = zeros(size(xx));
evalmat = chunkerkernevalmat(chnkr,ekern,targout,struct('corrections',1));
uscat(out) = pcfft_apply(sol, A_spread_s, A_spread_t,evalmat+A_addsub_eval,skern_hat);
tplot = toc

% uscat(out) = chunkerkerneval(chnkr,ekern,sol2,targout);
utot = uin + uscat;

%% make plots
umax = max(abs(utot(:))); 
figure(2);clf
h = pcolor(xx,yy,imag(utot)); set(h,'EdgeColor','none'); colorbar
colormap(redblue); clim([-umax,umax]);
hold on 
plot(chnkr,'k')
axis equal
title('$u^{\textrm{tot}}$','Interpreter','latex','FontSize',12)

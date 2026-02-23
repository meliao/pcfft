%Chunkie integration demo
%
% Solve Dirichlet scattering problems with many inclusions


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
theta = 2*pi*rand();
ctrs = L*rand()*[cos(theta);sin(theta)];
n_pts = [];
nscat = 25;
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

dkern = kernel('helm', 'd', zk);
skern = kernel('helm', 's', zk);
% spkern = kernel('helm', 'sp', zk);

% define boundary data
rhs = -planewave(kvec,chnkr.r(:,:));
rhs = rhs(:);
%% Compute example field

L = 1.1*max(vecnorm(chnkr.r(:,:)));
x1 = linspace(-L,L,500);
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

%% Precompuation

eps = 1e-6;

srcs = [];
srcs.r = chnkr.r(:,:);
srcs.n = chnkr.n(:,:);

t1 = tic;

% setup grid
[grid_info, proxy_info] = get_grid(skern.eval, chnkr, targout, eps);
% get spreading operators
[A_spread_s, ~, sort_info_c]= get_spread(skern.eval, dkern.eval, chnkr, ...
    grid_info, proxy_info, {'r','n'});
[A_spread_c, ~, ~]= get_spread(skern.eval, skern.eval, chnkr, ...
    grid_info, proxy_info);
[A_spread_t, ~, sort_info_t] = get_spread(skern.eval, skern.eval, targout, ...
    grid_info, proxy_info);
% build corrections
[A_addsub_c] = get_addsub(skern.eval, dkern.eval, chnkr, chnkr, ...
    grid_info, proxy_info, sort_info_c, sort_info_c, A_spread_s, A_spread_c);
[A_addsub_eval] = get_addsub(skern.eval, dkern.eval, chnkr, targout, ...
    grid_info, proxy_info, sort_info_c, sort_info_t, A_spread_s, A_spread_t);
% get DFT of kernel
skern_hat = get_kernhat(skern.eval,grid_info);

cors = chunkermat(chnkr,dkern,struct('corrections',1));
cors = cors + A_addsub_c.*chnkr.wts(:).' + 0.5*speye(size(cors));
A_spread_s = A_spread_s.*chnkr.wts(:).';
A_addsub_eval = A_addsub_eval.*chnkr.wts(:).';
tprecom = toc(t1)

sys_app = @(dens) pcfft_apply(dens,A_spread_s,A_spread_c,cors,skern_hat);
% return
%%
tic;
% build fast direct solver
F = chunkerflam(chnkr,dkern,0.5);
t_build_solver = toc

tic;
% solve
sol2 = rskelf_sv(F,rhs);
tsolve1 = toc

tic;
% solve
sol = gmres(sys_app,rhs,[],eps,1000);
tsolve2 = toc

corsfmm = chunkermat(chnkr,dkern,struct('corrections',1));
corsfmm = corsfmm +  0.5*speye(size(cors));
sys_app = @(dens) chunkermatapply(chnkr,dkern,dens,corsfmm);
tic;
% solve
sol = gmres(sys_app,rhs,[],eps,1000);
tsolve3 = toc

%%

% get solution
tic;
uscat = zeros(size(xx));
evalmat = chunkerkernevalmat(chnkr,dkern,targout,struct('corrections',1));
uscat(out) = pcfft_apply(sol, A_spread_s, A_spread_t,evalmat+A_addsub_eval,skern_hat);
tplot = toc

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

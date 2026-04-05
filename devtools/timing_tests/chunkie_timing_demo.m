% Chunkie timing demo
% Requires chunkIE: https://github.com/fastalgorithms/chunkie
%
% Solve Dirichlet scattering problems with many inclusions

% Need to load planewave function
addpath("utils")
warning('off','MATLAB:rankDeficientMatrix')
warning('off','MATLAB:nearlySingularMatrix')
eps = 1e-6;

% planewave direction
phi = 0;
% ambient wavenumber
zk = 10;
kvec = zk*[cos(phi);sin(phi)];

nscats = 5:5:30;
% nscats = 450;
times = zeros(3,length(nscats));
solvetimes = zeros(3,length(nscats));
npts = zeros(1,length(nscats));
for k = 1:length(nscats)

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
L = 9;L = 35;
theta = 2*pi*rand();
ctrs = L*rand()*[cos(theta);sin(theta)];
n_pts = [];
nscat = nscats(k);
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

t1 = tic;

% setup grid
[grid_info, proxy_info] = get_grid(skern, chnkr, targout, eps,200);
% get spreading operators
[A_spread_s, sort_info_c]= get_spread(skern, dkern, chnkr, ...
    grid_info, proxy_info, {'r','n'});
[A_spread_c, ~]= get_spread(skern, skern, chnkr, ...
    grid_info, proxy_info);
[A_spread_t, sort_info_t] = get_spread(skern, skern, targout, ...
    grid_info, proxy_info);
% build corrections
A_addsub_c = get_addsub(skern, dkern, grid_info, proxy_info, ...
    sort_info_c, sort_info_c, A_spread_s, A_spread_c);
A_addsub_eval = get_addsub(skern, dkern, grid_info, proxy_info, ...
    sort_info_c, sort_info_t, A_spread_s, A_spread_t);
% get DFT of kernel
skern_hat = get_kernhat(skern,grid_info);

cors = chunkermat(chnkr,dkern,struct('corrections',1));
cors = cors + A_addsub_c.*chnkr.wts(:).' + 0.5*speye(size(cors));
A_spread_s = A_spread_s.*chnkr.wts(:).';
A_addsub_eval = A_addsub_eval.*chnkr.wts(:).';
tpcfftprecom = toc(t1);

sys_app = @(dens) pcfft_apply(dens,A_spread_s,A_spread_c,cors,skern_hat);
% return
tic;
% solve
sol = gmres(sys_app,rhs,[],eps,1000);
tpcfftsolve = toc;

%%
tic;
% build fast direct solver
flam_opts = [];flam_opts.rank_or_tol = eps;
F = chunkerflam(chnkr,dkern,0.5,flam_opts);
tFLAMprecom = toc;

tic;
% solve
sol2 = rskelf_sv(F,rhs);
tflamsolve = toc;


tic;
corsfmm = chunkermat(chnkr,dkern,struct('corrections',1));
corsfmm = corsfmm +  0.5*speye(size(cors));
sys_app = @(dens) chunkermatapply(chnkr,dkern,dens,corsfmm);
tfmmprecom = toc;
tic;
% solve
sol = gmres(sys_app,rhs,[],eps,1000);
tfmmsolve = toc;

%%
fprintf('PCFFT precom took %.2e s. The solve took %.2e s.\n',tpcfftprecom,tpcfftsolve)
fprintf('FLAM precom took  %.2e s. The solve took %.2e s.\n',tFLAMprecom,tflamsolve)
fprintf('FMM precom took  %.2e s. The solve took %.2e s.\n',tfmmprecom,tfmmsolve)

fprintf('PCFFT total: %.2e s. FLAM total: %.2e s. FMM total: %.2e s.\n',tpcfftprecom+tpcfftsolve,tFLAMprecom+tflamsolve,tfmmsolve)

times(1,k) = tpcfftprecom+tpcfftsolve;
times(2,k) = tFLAMprecom+tflamsolve;
times(3,k) = tfmmprecom+tfmmsolve;

solvetimes(1,k) = tpcfftsolve;
solvetimes(2,k) = tflamsolve;
solvetimes(3,k) = tfmmsolve;
npts(k) = chnkr.npt;

end


%% compute utot on the 2D domain for plotting

% get solution
tic;
uscat = (NaN+NaN*1i)*zeros(size(xx));
evalmat = chunkerkernevalmat(chnkr,dkern,targout,struct('corrections',1));
uscat(out) = pcfft_apply(sol, A_spread_s, A_spread_t,evalmat+A_addsub_eval,skern_hat);
tplot = toc

utot = uin + uscat;

%% make plots
umax = max(abs(utot(:))); 
figure(2);clf
h = pcolor(xx,yy,imag(utot)); set(h,'EdgeColor','none'); colorbar
colormap(redblue); clim(0.4*[-umax,umax]);
hold on 
plot(chnkr,'k')
axis equal
% title('$u^{\textrm{tot}}$','Interpreter','latex','FontSize',12)


figure(4);
% plot(npts,times,'o-')
plot(log10(npts),log10(times),'o-','LineWidth',2)
hold on
plot(log10(npts),log10(npts/600),'k--','LineWidth',2)
plot(log10(npts),log10(npts/400).^(3/2),'k:','LineWidth',2)
hold off
xlabel('$\log_{10} n_{pts}$','interpreter','latex')
ylabel('$\log_{10} $ time (s)','interpreter','latex')
legend('PCFFT', 'FLAM', 'FMM', '$O(n)$', '$O(n^{3/2})$','interpreter','latex','Location','northwest')
set(gca,'ticklabelinterpreter','latex')
set(gca,'fontsize',16)
% exportgraphics(gcf,'chunkie_scatterer_timings_tol1e-6.pdf')

figure(5);
% plot(npts,times,'o-')
plot(log10(npts),log10(solvetimes),'o-','LineWidth',2)
hold on
plot(log10(npts),log10(npts/600),'k--','LineWidth',2)
plot(log10(npts),log10(npts/400).^(3/2),'k:','LineWidth',2)
hold off
xlabel('$\log_{10} n_{pts}$','interpreter','latex')
ylabel('$\log_{10} $ solve time (s)','interpreter','latex')
legend('PCFFT', 'FLAM', 'FMM', '$O(n)$', '$O(n^{3/2})$','interpreter','latex','Location','northwest')
set(gca,'ticklabelinterpreter','latex')
set(gca,'fontsize',16)
% exportgraphics(gcf,'chunkie_scatterer_solvetimings_tol1e-6.pdf')
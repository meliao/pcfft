

nover = 2; L = 40;
nch = ceil(36*2^(nover-1));

cparams = [];
cparams.ta = -L-8;
cparams.tb =  L+8;

% cparams.ta = -L+5;
% cparams.tb =  L-5;

cpars = [];
cpars.L = L;
cpars.c2 = 5;
cpars.c1 = 20;
chnkr = chunkerfuncuni(@(t) complexx2(t,cpars),...
   nch,cparams);

alpha = 1;
beta = 100;

skern = kernel('laplace','s');

% kkern_fun = @(s,t) exp(1i*sqrt(beta/alpha) * sqrt((s.r(1,:)-t.r(1,:).').^2)) ...
%     * sqrt(alpha/beta)/2i;

kkern_fun = @(s,t) chnk.helm1d.kern(sqrt(beta/alpha), s,t,'s');
kkern = kernel(kkern_fun);


tic;
Smat = chunkermat(chnkr,skern);
Kmat = chunkermat(chnkr,kkern);
toc;

tic;
sysmat = eye(chnkr.npt) + (2/alpha)*Smat*Kmat;
toc;

src = [0;0];
f = exp(-vecnorm(chnkr.r(:,:)-src).^2/10).';

tic;
mu = sysmat\(2/alpha * f);
toc;

sigma = Kmat*mu;

%%
figure(1);clf
subplot(1,3,1)
plot(real(chnkr.r(1,:)), real(mu),'.')
hold on,
plot(real(chnkr.r(1,:)), imag(mu),'.')
hold off

subplot(1,3,2)
plot(real(chnkr.r(1,:)), log10(abs(mu)),'.')

subplot(1,3,3)
plot(real(chnkr.r(1,:)), real(sigma),'.')
hold on,
plot(real(chnkr.r(1,:)), imag(sigma),'.')
hold off


% figure(2)
% plot(real(chnkr.r(1,:)),f)









% Compare dx_nproxy against dx_nproxy_old across kernels and tolerances.
%
% Reports the parameters each one returns, the time taken, and the error those
% parameters actually achieve on a fresh random configuration.

addpath(genpath('../../pcfft'));

halfside = 0.5;
crad     = 2;
nrep     = 3;          % repeats; the minimum time is reported

cases = { ...
    struct('name','log 2D',        'kern',@log_kernel,          'dim',2), ...
    struct('name','one_over_r 2D', 'kern',@one_over_r_kernel2D, 'dim',2), ...
    struct('name','log 3D',        'kern',@log_kernel3D,        'dim',3), ...
    struct('name','one_over_r 3D', 'kern',@one_over_r_kernel,   'dim',3)};

tols = 10.^(-(3:2:13));

for ic = 1:numel(cases)
    c = cases{ic};
    if isempty(which(func2str(c.kern)))
        fprintf('\n=== %s (not on path, skipped) ===\n', c.name);
        continue
    end

    fprintf('\n=== %s ===\n', c.name);
    fprintf('%8s | %-31s | %-31s | %8s\n', ...
            'tol', 'NEW  nspread nproxy nshell     t', ...
                   'OLD  nspread nproxy nshell     t', 'speedup');
    fprintf('%s\n', repmat('-',1,88));

    for tol = tols
        [new, msg_new] = run_one(@dx_nproxy,     c, tol, halfside, crad, nrep);
        [old, msg_old] = run_one(@dx_nproxy_old, c, tol, halfside, crad, nrep);

        fprintf('%8.1e | %7d %6d %6d %7.3fs | %7d %6d %6d %7.3fs | %8s\n', ...
                tol, new.nspread, new.nproxy, new.nshell, new.t, ...
                     old.nspread, old.nproxy, old.nshell, old.t, ...
                     speedup(old.t, new.t));
        if ~isempty(msg_new), fprintf('%10s NEW FAILED: %s\n', '', msg_new); end
        if ~isempty(msg_old), fprintf('%10s OLD FAILED: %s\n', '', msg_old); end
        if ~isnan(new.err) || ~isnan(old.err)
            fprintf('%10s observed err  new %8.2e   old %8.2e\n', '', new.err, old.err);
        end
    end
end

fprintf('\nDone.\n');


%% ------------------------------------------------------------------
function [r, msg] = run_one(fn, c, tol, halfside, crad, nrep)
    r = struct('nspread',NaN,'nproxy',NaN,'nshell',NaN,'t',NaN,'err',NaN);
    msg = '';
    try
        t = Inf;
        for k = 1:nrep
            t0 = tic;
            [~, nspread, ~, pinfo] = fn(c.kern, c.dim, tol, halfside, crad);
            t = min(t, toc(t0));
        end
        r.t = t;
        r.nspread = nspread;
        r.nproxy  = pinfo.nproxy;
        r.nshell  = pinfo.nshell;
        r.err     = check_error(c.kern, c.dim, halfside, nspread, pinfo);
    catch ME
        msg = ME.message;
    end
end

function s = speedup(to, tn)
    if isnan(to) || isnan(tn) || tn == 0
        s = '-';
    else
        s = sprintf('%.1fx', to/tn);
    end
end

function err = check_error(kern, dim, halfside, nspread, pinfo)
    % Error achieved by these parameters on a fresh random configuration.
    rng(12345);
    radius = pinfo.radius;
    bwidth = halfside;

    nsrc = 200;
    src  = (rand(dim,nsrc)-0.5)*bwidth;
    w    = rand(nsrc,1)-0.5;

    d = randn(dim,300); d = d./vecnorm(d);
    targ = d .* (radius*(1.01 + 3*rand(1,300)));

    K = @(s,t) kern(struct('r',s), struct('r',t));

    dx = 2*halfside/nspread;
    xx = -halfside + dx/2 + (0:nspread-1)*dx;
    if dim == 2
        [X,Y] = meshgrid(xx,xx);
        reg = [X(:).';Y(:).'];
    else
        [X,Y,Z] = meshgrid(xx,xx,xx);
        X=permute(X,[3,1,2]); Y=permute(Y,[3,1,2]); Z=permute(Z,[3,1,2]);
        reg = [X(:).';Y(:).';Z(:).'];
    end

    Kp = @(s,t) wrap_kern_der(kern, struct('r',s), struct('r',t), pinfo.proxy_der);
    sw = lsqminnorm(Kp(reg, pinfo.r), Kp(src, pinfo.r)*w, pinfo.tol/10);

    exact = K(src,targ)*w;
    err   = max(abs(K(reg,targ)*sw - exact))/max(abs(exact));
end

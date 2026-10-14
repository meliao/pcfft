function [dx, nspread, nbinpts, proxy_info] = dx_nproxy(kernel, dim, tol, halfside, crad, multi_shells, proxy_der, opts)
% DX_NPROXY Determines the grid spacing dx, number of spreading points
% across the spreading box nspread, number of grid points across the
% spreading bin nbinpts, and proxy surface information proxy_info. We
% do this in the following steps:
%
%   1. Find nproxy that resolves K(.,x) down to tol on a single shell
%   by looking at Fourier or spherical harmonic tails. Leaves a 25%
%   margin in the tails.
%
%   2. For increasing grid sizes, try to find allowable nshell. We take the
%   first (nspread, nshell) pair that achieves the tolerance at targets.
%
% Parameters
% ----------
% kernel : function handle or object with .eval
% dim : dimension
% tol : Relative L_inf error tolerance.
% halfside : Halfside of the spreading box.
% crad : radius of the innermost proxy shell is = crad * sqrt(dim) * halfside.
% multi_shells : allow more than one proxy shell. Defaults to true; false
%   locks nshell at 1.
% proxy_der : Number of derivatives used in proxy compression. 
% opts : struct of search limits, any subset of
%     opts.nproxy_max
%         Maximum number of proxy points per shell (default 6000)
%     opts.nshell_max
%         Maximum number of shells (default 20)
%     opts.nspread_max
%         Maximum grid points across the box (default 100 in 2D, 30 in 3D)
%     opts.relax_nshell
%         Trade one more nspread for a smaller nshell (default false)
%     opts.spread_cases
%         Cell array of structs describing additional (derivative) kernels
%         that will be spread with this grid, each with fields
%           .kern   : kernel, called as kern(pt, struct('r', x)) where pt
%                     carries the derivative fields (same convention as
%                     the kern_der argument of get_spread)
%           .info   : point struct the fields are sampled from
%           .fields : cell array of field names to sample from .info
%         nproxy is taken as the max over all cases, and a grid is only
%         accepted if every case meets the tolerance. The free-space
%         kernel is always included. Usually filled in by get_grid.
%
% Returns
% -------
% dx, nspread, nbinpts, proxy_info

    if nargin < 5 || isempty(crad)
        crad = 2;
    end
    if nargin < 6 || isempty(multi_shells)
        multi_shells = true;
    end
    if nargin < 7 || isempty(proxy_der)
        proxy_der = 0;
    end
    if nargin < 8
        opts = struct();
    end

    if ~isa(kernel, 'function_handle')
        try
            kernel = kernel.eval;
        catch
            error('dx_nproxy: kernel is not a function and does not have an eval property')
        end
    end

    % Grid -> proxy / target always uses the free-space kernel.
    kernel_pxy = @(s, t) wrap_kern_der(kernel, struct('r', s), struct('r', t), proxy_der);
    kernel_ev  = @(s, t) kernel(struct('r', s), struct('r', t));

    % Kernels that will actually be spread. The free-space kernel is always
    % one of them; derivative kernels are added from opts.spread_cases.
    cases = {struct('kern', kernel, 'info', [], 'fields', {{}})};
    if isfield(opts, 'spread_cases') && ~isempty(opts.spread_cases)
        extra = opts.spread_cases;
        if isstruct(extra), extra = num2cell(extra); end
        for k = 1:numel(extra)
            c = extra{k};
            if ~isa(c.kern, 'function_handle')
                try
                    c.kern = c.kern.eval;
                catch
                    error('dx_nproxy: spread_cases{%d}.kern is not a function and does not have an eval property', k)
                end
            end
            if ~isfield(c, 'fields') || isempty(c.fields), c.fields = {}; end
            if ischar(c.fields), c.fields = {c.fields}; end
            if ~isfield(c, 'info'), c.info = []; end
            if ~isempty(setdiff(c.fields, {'r'})) && isempty(c.info)
                error('dx_nproxy: spread_cases{%d} requests fields but has no info to sample them from', k)
            end
            cases{end+1} = c; %#ok<AGROW>
        end
    end
    ncase = numel(cases);

    rng(0);

    % Geometry
    c_bwidth = 1.0;
    bwidth   = c_bwidth * halfside;

    R      = sqrt(dim) * halfside;
    radius = crad * R;

    % The solve costs nspread^dim, so 3D gives up sooner.
    if dim == 2
        nspread_max_default = 100;
    else
        nspread_max_default = 30;
    end

    nproxy_max = 6000;
    if isfield(opts,'nproxy_max')
        nproxy_max = opts.nproxy_max;
    end
    nshell_max = 20;
    if isfield(opts,'nshell_max')
        nshell_max = opts.nshell_max;
    end
    if ~multi_shells
        nshell_max = 1;
    end
    nspread_max = nspread_max_default;
    if isfield(opts,'nspread_max')
        nspread_max = opts.nspread_max;
    end

    % Sources and targets for the acceptance tests.
    nsrc = 117;
    src_pts     = (rand(dim, nsrc) - 0.5) * bwidth;
    src_weights = rand(nsrc, 1) - 0.5;   % used by the free-space case

    % Spread the targets over an annulus
    ntarget = 500;
    if dim == 2
        target_pts = (1 + rand(1, ntarget)) .* get_ring_points(ntarget, 1.01 * radius);
    else
        target_pts = (1 + rand(1, ntarget)) .* get_sphere_points(ntarget, 1.01 * radius);
    end

    % Corners are added when checking the resolution of the field
    if dim == 2
        [A, B] = ndgrid([-1 1], [-1 1]);
        c = [A(:).'; B(:).'];
    else
        [A, B, C] = ndgrid([-1 1], [-1 1], [-1 1]);
        c = [A(:).'; B(:).'; C(:).'];
    end
    src_res_pts = [src_pts, c * (bwidth / 2)];
    nres = size(src_res_pts, 2);

    % Per-case test data. Fields are sampled once (after rng(0)) so every
    % evaluation in the search sees the same synthetic sources.
    tests = cell(1, ncase);
    for k = 1:ncase
        ck = cases{k};
        if isempty(ck.info)
            idx = ones(1, nres);
        else
            idx = randi(size(ck.info.r(:,:), 2), 1, nres);
        end
        S_res = sample_point_fields(src_res_pts, ck.info, ck.fields, idx);
        S     = sample_point_fields(src_pts, ck.info, ck.fields, idx(1:nsrc));

        T.src_pxy = @(t) wrap_kern_der(ck.kern, S, struct('r', t), proxy_der);
        T.res_pxy = @(t) wrap_kern_der(ck.kern, S_res, struct('r', t), proxy_der);
        K_st = ck.kern(S, struct('r', target_pts));
        opdim = size(K_st, 2) / nsrc;
        if k == 1 && opdim == 1
            T.weights = src_weights;
        else
            T.weights = rand(opdim * nsrc, 1) - 0.5;
        end
        T.target_evals = K_st * T.weights;
        T.scale_targ = max(abs(T.target_evals));
        tests{k} = T;
    end

    % Stage 1: Find nproxy that resolves the field of every kernel
    nproxy = 0;
    for k = 1:ncase
        if dim == 2
            npk = nproxy_ring(tests{k}.res_pxy, radius, tol, nproxy_max);
        else
            npk = nproxy_sphere(tests{k}.res_pxy, radius, tol, nproxy_max);
        end
        nproxy = max(nproxy, npk);
    end

    % Stages 2 and 3: Find nshell and nspread
    % For each grid size, increase nshell looking for allowable parameters

    relax_nshell = false;
    if isfield(opts,'relax_nshell')
        relax_nshell = opts.relax_nshell;
    end

    if dim == 2
        nmin = 4;
    else
        nmin = 3;
    end

    nshell  = [];
    nspread = [];
    for n = nmin:nspread_max
        for ns = 1:nshell_max
            proxy_pts = proxy_surface(dim, nproxy, ns, radius);

            err = spread_error(kernel_pxy, kernel_ev, tests, dim, n, halfside, ...
                               proxy_pts, target_pts, tol);

            if err < tol
                nshell  = ns;
                nspread = n;
                dx      = 2 * halfside / n;
                break
            end
        end

        if ~isempty(nshell)
            break
        end

        % No shell count worked at this grid size: refine the grid.
    end
    if isempty(nshell)
        if ~multi_shells
            error('dx_nproxy: no grid fond for one shell, try multi_shells = true')
        end
        error('dx_nproxy: no grid up to nspread = %d works with up to %d shells', ...
              nspread_max, nshell_max)
    end

    % Test if one more grid point yields fewer shells.
    if relax_nshell && nshell > 1 && nspread < nspread_max
        n = nspread + 1;
        for ns = 1:nshell - 1
            proxy_pts = proxy_surface(dim, nproxy, ns, radius);

            err = spread_error(kernel_pxy, kernel_ev, tests, dim, n, halfside, ...
                               proxy_pts, target_pts, tol);

            if err < tol && ns < nshell
                nshell  = ns;
                nspread = n;
                dx      = 2 * halfside / n;
                break
            end
        end
    end

    proxy_pts = proxy_surface(dim, nproxy, nshell, radius);

    % Assemble the output, avoiding a square least squares matrix
    if nproxy * nshell == nspread^dim
        nproxy    = nproxy + 1;
        proxy_pts = proxy_surface(dim, nproxy, nshell, radius);
    end

    nbinpts = floor(nspread / 2);

    proxy_info = ProxyInfo(dim, nproxy * nshell, nproxy, nshell, halfside, ...
                           crad, tol, radius, proxy_pts, proxy_der);

end

function err = spread_error(kernel_pxy, kernel_ev, tests, dim, nspread, halfside, ...
                            proxy_pts, target_pts, tol)
    % Measure the worst approximation error over all spread kernels
    % obtained with this grid
    dx = 2 * halfside / nspread;
    xx = -halfside + dx / 2 + (0:nspread - 1) * dx;
    if dim == 2
        [X, Y] = meshgrid(xx, xx);
        reg_pts = [X(:).'; Y(:).'];
    else
        [X, Y, Z] = meshgrid(xx, xx, xx);
        X = permute(X, [3, 1, 2]);
        Y = permute(Y, [3, 1, 2]);
        Z = permute(Z, [3, 1, 2]);
        reg_pts = [X(:).'; Y(:).'; Z(:).'];
    end

    K_reg_to_proxy = kernel_pxy(reg_pts, proxy_pts);
    K_reg_to_targ  = kernel_ev(reg_pts, target_pts);

    err = 0;
    for k = 1:numel(tests)
        T = tests{k};
        evals_at_proxy = T.src_pxy(proxy_pts) * T.weights;
        spread_weights = lsqminnorm(K_reg_to_proxy, evals_at_proxy, tol / 10);
        approx_at_target = K_reg_to_targ * spread_weights;
        errk = max(abs(approx_at_target(:) - T.target_evals(:))) / T.scale_targ;
        err = max(err, errk);
        if err >= tol, return; end
    end
end

function nproxy = nproxy_ring(field_at, radius, tol, nproxy_max)
    % field_at(p) returns the (derivative-wrapped) kernel matrix from the
    % test sources to the points p.
    % Find point where Fourier tails of k(.,src) fall below tol. Leaves a
    % margin of at least 3 orders and at least 25%
    nmarg_abs  = 3;
    nmarge_frac = 0.25;

    % Resolve every mode the sweep can reach.
    nfine = 2 * nproxy_max;
    p     = get_ring_points(nfine, radius);
    F     = field_at(p);

    nder  = size(F, 1) / nfine;
    nh    = floor(nfine / 2);
    tails = zeros(nh, 1);
    % treat each derivative as its own function of the angle.
    for k = 1:nder
        A = abs(fft(F(k:nder:end, :), [], 1) / nfine);

        % fold the negative frequencies onto the positive ones
        B = A(1:nh, :);
        B(2:nh, :) = max(B(2:nh, :), A(nfine:-1:nfine - nh + 2, :));

        % running max from each mode outward, each source scaled by its peak
        T = flipud(cummax(flipud(B), 1)) ./ max(B, [], 1);

        tails = max(tails, max(T, [], 2));
    end

    for np = 8:4:nproxy_max
        mnyq = floor(np / 2);
        m    = mnyq - max(nmarg_abs, ceil(nmarge_frac * mnyq));
        if m < 1
            continue
        end
        if tails(m + 1) < tol
            nproxy = np;
            return
        end
    end
    error('dx_nproxy: nproxy did not converge below nproxy = %d', nproxy_max)
end

function nproxy = nproxy_sphere(field_at, radius, tol, nproxy_max)
    % field_at(p) returns the (derivative-wrapped) kernel matrix from the
    % test sources to the points p.
    % Find point where spherical harmonic tails of k(.,src) fall below tol. Leaves a
    % margin of at least 3 orders and at least 25%
    nmarg_abs  = 3;
    nmarge_frac = 0.25;
    n_for_degree = @(L) 1.05 * (L + 1)^2;

    % Resolve every degree the sweep can reach.
    Lmax   = ceil(sqrt(nproxy_max / 1.05));
    ntheta = Lmax + 1;
    nphi   = 2 * Lmax + 2;

    [ct, wt] = gauss_legendre(ntheta);
    st       = sqrt(1 - ct(:).^2);
    phi      = (0:nphi - 1) * (2 * pi / nphi);

    X = radius * (st * cos(phi));
    Y = radius * (st * sin(phi));
    Z = radius * repmat(ct(:), 1, nphi);
    p = [X(:).'; Y(:).'; Z(:).'];

    ngrid = ntheta * nphi;
    F     = field_at(p);
    nder  = size(F, 1) / ngrid;
    nsrc  = size(F, 2);

    % Per-degree energy, kept per source so each is normalized by its own peak.
    deg_energy = zeros(Lmax + 1, nsrc);
    for k = 1:nder
        Fk = reshape(F(k:nder:end, :), ntheta, nphi, nsrc);
        G  = fft(Fk, [], 2) / nphi;

        for l = 0:Lmax
            P = legendre(l, ct(:).', 'norm');
            e = zeros(1, nsrc);
            for m = 0:l
                w  = (wt(:) .* P(m + 1, :).');
                cm = reshape(sum(G(:, m + 1, :) .* w, 1), 1, nsrc);
                e  = max(e, abs(cm));
                if m > 0
                    cn = reshape(sum(G(:, nphi - m + 1, :) .* w, 1), 1, nsrc);
                    e  = max(e, abs(cn));
                end
            end
            deg_energy(l + 1, :) = max(deg_energy(l + 1, :), e);
        end
    end

    T    = flipud(cummax(flipud(deg_energy), 1));
    T    = T ./ max(deg_energy, [], 1);
    tail = max(T, [], 2);

    for L = 0:Lmax - 1
        Lc = L - max(nmarg_abs, ceil(nmarge_frac * L));
        if Lc < 0
            continue
        end
        if tail(Lc + 2) < tol
            nproxy = ceil(n_for_degree(L));
            if nproxy > nproxy_max
                break
            end
            return
        end
    end
    error('dx_nproxy: nproxy did not converge below nproxy = %d', nproxy_max)
end

function [x, w] = gauss_legendre(n)
    % Golub-Welsch
    k = 1:n - 1;
    b = k ./ sqrt(4 * k.^2 - 1);
    J = diag(b, 1) + diag(b, -1);
    [V, D] = eig(J);
    [x, i] = sort(diag(D));
    w = 2 * (V(1, i).^2).';
    x = x(:);
    w = w(:);
end


function proxy_pts = proxy_surface(dim, nproxy, nshell, radius)
    if dim == 2
        proxy_pts0 = get_ring_points(nproxy, 1);
    else
        proxy_pts0 = get_sphere_points(nproxy, 1);
    end

    % Chebyshev-like spacing in 1/r, starting from the innermost shell.
    if nshell > 1
        sr   = cos((1:nshell - 1) / nshell * pi);
        sr   = (sr(:) + 1) / 2 * (1 / radius);
        rads = [radius, 1 ./ sr(:).'];
    else
        rads = radius;
    end

    proxy_pts = reshape(rads(:).' .* reshape(proxy_pts0, dim, 1, []), dim, []);
end

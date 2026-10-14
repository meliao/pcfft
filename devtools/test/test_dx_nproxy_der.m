addpath(genpath('../../pcfft'));

% Source-derivative kernels need more proxy points than the free-space
% kernel. Check that passing them through opts increases nproxy and
% decreases error.
%
% Test kernels are derivatives of log|t - s| in the source normal
% direction nu.

tol = 1e-8;
rng(1);

kern_0 = @log_kernel;
cz  = @(s, t) (t.r(1,:).' - s.r(1,:)) + 1i*(t.r(2,:).' - s.r(2,:));
cnu = @(s) s.n(1,:) + 1i*s.n(2,:);
dlp = @(s, t) real(-cnu(s)    ./ cz(s, t));      % 1st normal derivative
d3  = @(s, t) real(-2*cnu(s).^3 ./ cz(s, t).^3); % 3rd normal derivative

kerns = {dlp, d3};
names = {'dlp', 'd3'};

% Sources on a circle carry normals; targets carry nothing extra.
nsrc = 2000;
th = 2*pi*rand(1, nsrc);
src_info  = struct('r', [cos(th); sin(th)], 'n', [cos(th); sin(th)]);
targ_info = struct('r', 2*rand(2, 1500) - 1);

[~, proxy_0] = get_grid(kern_0, src_info, targ_info, tol);
nproxy_prev = proxy_0.nproxy;

far = struct('r', 6 * [cos(2*pi*(0:99)/100); sin(2*pi*(0:99)/100)]);
mu  = rand(nsrc, 1) - 0.5;

for k = 1:numel(kerns)
    kd = kerns{k};

    %% nproxy grows with the order of the derivative
    opts = struct('kern_der_src', kd, 'der_fields_src', {{'r', 'n'}});
    [grid_info, proxy_d] = get_grid(kern_0, src_info, targ_info, tol, [], opts);
    fprintf('%s: nproxy = %d (free space %d), nshell = %d, nspread = %d\n', ...
        names{k}, proxy_d.nproxy, proxy_0.nproxy, proxy_d.nshell, grid_info.nspread);
    assert(proxy_d.nproxy > nproxy_prev, ...
        '%s: expected nproxy > %d, got %d', names{k}, nproxy_prev, proxy_d.nproxy);
    nproxy_prev = proxy_d.nproxy;

    %% default fields are inferred from src_info
    [~, proxy_inf] = get_grid(kern_0, src_info, targ_info, tol, [], ...
                              struct('kern_der_src', kd));
    assert(proxy_inf.nproxy == proxy_d.nproxy);

    %% the grid spreads the derivative kernel to tolerance at far targets
    [A_spread, sort_info] = get_spread(kern_0, kd, src_info, grid_info, ...
                                       proxy_d, {'r', 'n'});
    u_ref = kd(src_info, far) * mu;
    u_pc  = kern_0(struct('r', grid_info.r), far) * (A_spread * mu(sort_info.ptid_srt));
    err = max(abs(u_pc - u_ref)) / max(abs(u_ref));
    fprintf('%s: far-field spread error %.2e\n', names{k}, err);
    assert(err < 100*tol, '%s spread error %.2e', names{k}, err);

    %% without the derivative kernel in opts, the grid is under-resolved
    [grid_0, proxy_00] = get_grid(kern_0, src_info, targ_info, tol);
    [A0, si0] = get_spread(kern_0, kd, src_info, grid_0, proxy_00, {'r', 'n'});
    u0 = kern_0(struct('r', grid_0.r), far) * (A0 * mu(si0.ptid_srt));
    err0 = max(abs(u0 - u_ref)) / max(abs(u_ref));
    fprintf('%s: error with free-space-only grid %.2e\n', names{k}, err0);
end

%% fields for target kernels come from targ_info, so asking for 'n' fails
opts_bad = struct('kern_der_targ', d3, 'der_fields_targ', {{'r', 'n'}});
threw = false;
try
    get_grid(kern_0, src_info, targ_info, tol, [], opts_bad);
catch
    threw = true;
end
assert(threw, 'expected an error for a target field that does not exist');

%% target-side kernel with target normals
targ_n = targ_info;
phi = 2*pi*rand(1, size(targ_n.r, 2));
targ_n.n = [cos(phi); sin(phi)];
opts_t = struct('kern_der_targ', d3, 'der_fields_targ', {{'r', 'n'}});
[~, proxy_t] = get_grid(kern_0, src_info, targ_n, tol, [], opts_t);
assert(proxy_t.nproxy > proxy_0.nproxy);

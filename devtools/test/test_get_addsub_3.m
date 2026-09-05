% Exercises get_addsub in the regime where the TARGETS occupy many more bins
% than the SOURCES: the sources sit in one tiny cluster while the targets fill
% a much larger box. Both point sets share a single bin decomposition built by
% get_grid, so the bin loop in get_addsub -- whose bound comes from
% sort_info_s.id_start but which indexes sort_info_t.id_start(i) -- must still
% visit every occupied target bin.
%
% Correctness is checked against the dense equivalent of pcfft_apply,
%
%     kern_0(src, targ)  ~=  A_addsub + A_spread_t.' * K_gg * A_spread_s,
%
% where K_gg is the grid-to-grid kernel with the self-interaction (diagonal)
% zeroed, matching the kernvals(1,1) = 0 in get_kernhat. The identity is exact
% for bin pairs inside the neighbour stencil and accurate to the spreading
% tolerance for the rest, so the assertions compare against tol.
%
% NOTE: sources and targets must stay disjoint (get_addsub zeroes coincident
% pairs but log_kernel does not), and the point coordinates must be random --
% a point sitting exactly on the upper face of the bounding box would be sorted
% into bin index nbin by SortInfo.
addpath(genpath('../../pcfft'));
warning('off', 'MATLAB:rankDeficientMatrix');

%% Source cluster in an interior bin
rng(0);
kern_0 = @(s,t) log_kernel(s,t);

tol  = 1e-6;
opts = struct('halfside', 0.13);   % force many bins; the bin side is halfside here

targ_info = struct('r', 2*rand(2, 1500) - 1);                % fills [-1, 1]^2
src_info  = struct('r', 0.01*randn(2, 200) + [0.13; -0.07]); % tiny interior cluster

[grid_info, proxy_info] = get_grid(kern_0, src_info, targ_info, tol, [], opts);

% Guard the dense reference: K_gg is nreg x nreg.
nreg = prod(grid_info.ngrid);
assert(nreg < 4000, 'grid too large for a dense reference (nreg = %d)', nreg);

[A_spread_s, sort_info_s, spread_blk_s] = get_spread(kern_0, [], src_info, ...
    grid_info, proxy_info);
[A_spread_t, sort_info_t, spread_blk_t] = get_spread(kern_0, [], targ_info, ...
    grid_info, proxy_info);

% Premise of the test. Both SortInfos index the same bins (get_addsub relies on
% this), but the targets occupy far more of them than the sources do.
assert(numel(sort_info_s.id_start) == numel(sort_info_t.id_start), ...
    'source and target SortInfo must share the same bin count');
n_occ_s = nnz(diff(sort_info_s.id_start));
n_occ_t = nnz(diff(sort_info_t.id_start));
assert(n_occ_t > 20 * n_occ_s, ...
    'premise: occupied target bins (%d) must far exceed occupied source bins (%d)', ...
    n_occ_t, n_occ_s);

A_addsub = get_addsub(kern_0, [], grid_info, proxy_info, sort_info_s, ...
    sort_info_t, spread_blk_s, spread_blk_t);

assert(all(size(A_addsub) == [size(targ_info.r,2), size(src_info.r,2)]));

% Dense reference: full grid-to-grid kernel with the diagonal zeroed.
K_gg = kern_0(struct('r', grid_info.r), struct('r', grid_info.r));
K_gg(1:nreg+1:end) = 0;    % assign, do not multiply: log(0) = -Inf and 0*-Inf = NaN

K_exact = kern_0(src_info, targ_info);
E   = full(A_addsub) + A_spread_t.' * (K_gg * A_spread_s) - K_exact;
err = max(abs(E(:))) / max(abs(K_exact(:)));
assert(err < tol, '2D interior cluster: rel err %g exceeds tol %g', err, tol);

%% Source cluster in the corner bin
% Same imbalance, but now the sources land in bin 0, so most of that bin's
% neighbour offsets are out of range and get filtered out of the correction.
rng(3);
kern_0 = @(s,t) log_kernel(s,t);

tol  = 1e-8;
opts = struct('halfside', 0.13);

targ_info = struct('r', 2*rand(2, 1500) - 1);
src_info  = struct('r', 0.008*randn(2, 150) + [-0.93; -0.90]);

[grid_info, proxy_info] = get_grid(kern_0, src_info, targ_info, tol, [], opts);

nreg = prod(grid_info.ngrid);
assert(nreg < 4000, 'grid too large for a dense reference (nreg = %d)', nreg);

[A_spread_s, sort_info_s, spread_blk_s] = get_spread(kern_0, [], src_info, ...
    grid_info, proxy_info);
[A_spread_t, sort_info_t, spread_blk_t] = get_spread(kern_0, [], targ_info, ...
    grid_info, proxy_info);

n_occ_s = nnz(diff(sort_info_s.id_start));
n_occ_t = nnz(diff(sort_info_t.id_start));
assert(n_occ_t > 20 * n_occ_s, ...
    'premise: occupied target bins (%d) must far exceed occupied source bins (%d)', ...
    n_occ_t, n_occ_s);
assert(sort_info_s.binid_srt(1) == 0, 'sources should land in the corner bin');

A_addsub = get_addsub(kern_0, [], grid_info, proxy_info, sort_info_s, ...
    sort_info_t, spread_blk_s, spread_blk_t);

K_gg = kern_0(struct('r', grid_info.r), struct('r', grid_info.r));
K_gg(1:nreg+1:end) = 0;

K_exact = kern_0(src_info, targ_info);
E   = full(A_addsub) + A_spread_t.' * (K_gg * A_spread_s) - K_exact;
err = max(abs(E(:))) / max(abs(K_exact(:)));
assert(err < tol, '2D corner cluster: rel err %g exceeds tol %g', err, tol);

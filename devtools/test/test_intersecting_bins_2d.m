% Makes sure get_addsub returns without error on a 2D input.
addpath(genpath('../../pcfft'));


rad = 10.0;
tol = 1e-10;
dim = 2;

k = @(s,t) log_kernel(s,t);

src_info_2d = struct;
n_src = 13;
rng(0);
src_info_2d.r = (rand(2, n_src) - 0.5);
src_info_2d.weights = rand(n_src, 1);

targ_info_2d = struct;
ntarg = 17;
targ_info_2d.r = rand(2, ntarg) - 0.5;


[grid_info, proxy_info] = get_grid(@log_kernel, ...
    src_info_2d, targ_info_2d, tol);



N_bin = grid_info.nbin(1) * grid_info.nbin(2);

% disp("test_intersecting_bins_2d: N_bin = " + int2str(N_bin));

% valid bin_idxes should be between 0 and N_bin - 1
for bin_idx = 0:(N_bin - 1)
    [idx_x, idx_y, binids] = intersecting_bins_2d(bin_idx, grid_info, proxy_info);

    % Once we remove invalid binids, there should be no repeats.
    valid_binids = binids(binids >= 0);
    assert(length(valid_binids) == length(unique(valid_binids)));

end

%% test_0b
% Larger test case with visualization. Same setup as test_SortInfo_2d_1
% TODO: Re-write this test once the nearness radius is figured out.

% n_pts = 10000;
% L = 2.0;
% Lbd = [-1 1; 
%         -1 1];
% % r points live on [-1, 1] x [-1, 1]
% rng(0);
% r = (rand(2, n_pts) - 0.5) * L;
% src_info_0b = struct('r', r);

% % Get grid and proxy info. The halfside is tuned so that the proxy shells of 
% % opposing corners do not intersect, see the figure generated.
% tol = 1e-6;
% [grid_info, proxy_info] = get_grid(@log_kernel, ...
%     src_info_0b, ...
%     src_info_0b, ...
%     tol, 1000, struct('halfside', 0.41));

% sort_info = SortInfo(src_info_0b, grid_info.dx, grid_info.Lbd, grid_info.nbin, grid_info.nbinpts);
% N_bins = grid_info.nbin(1) * grid_info.nbin(2);

% % Plot the sorted points and color by the bin
% % to make sure the bin assignment looks correct
% scatter(sort_info.r_srt(1,:), sort_info.r_srt(2,:), 20, sort_info.binid_srt, 'filled');
% colormap('parula');
% colorbar;

% % Draw an x at the center of each bin
% hold on;

% bin_penultimate = (grid_info.nbin(1) - 1) * grid_info.nbin(2) - 2;

% % Draw proxy rings for the first and last bin.
% for bin_idx =  [0 bin_penultimate N_bins - 1]
%     center = bin_center(bin_idx, grid_info);
%     scatter(center(1), center(2), 100, 'x');

%     % Draw the proxy ring
%     proxypts = get_ring_points(100, proxy_info.radius, center);
%     plot(proxypts(1,:), proxypts(2,:), 'k-');
% end


% % close all;


% disp("test_intersecting_bins_2d: grid_info:");
% disp(grid_info);
% disp("test_intersecting_bins_2d: grid_info.nbin:");
% disp(grid_info.nbin);

% % Test the third return value is correct.

% [~, ~, bin_0_intersecting_binids] = ...
%     intersecting_bins_2d(0, grid_info, proxy_info);

% % This should = [0 1 2 3 4 5 6 7 8 ... Nbin-1]
% expected_bin_0_intersecting_binids = 0:(N_bins - 1);
% % disp("test_intersecting_bins_2d: For bin_idx 4, intersecting binids: ");
% % disp(bin_4_intersecting_binids);
% valid_bins = bin_0_intersecting_binids >= 0 & bin_0_intersecting_binids < N_bins;
% valid_bins = bin_0_intersecting_binids(valid_bins);
% % disp("test_intersecting_bins_2d: For bin_idx 4, valid intersecting binids: ");
% % disp(valid_bins);
% assert(all(valid_bins == expected_bin_0_intersecting_binids));


% [sort_info] = SortInfo(struct('r', r), grid_info.dx, grid_info.Lbd, grid_info.nbin, grid_info.nbinpts);
% r_srt = sort_info.r_srt;
% binid_srt = sort_info.binid_srt;
% ptid_srt = sort_info.ptid_srt;
% id_start = sort_info.id_start;
% c = 1:n_pts;
% assert(all(size(r_srt) == size(r)));
% assert(all(size(r, 2) == size(binid_srt, 2)));

% % Figure shows that bin idx 0 intersercts with all of the bins in the first col.
% [bin_0_intersecting_x, bin_0_intersecting_y, ~] = intersecting_bins_2d(0, grid_info, proxy_info);

% % Expect bin_0_intersecting_x = [0, 1, ..., Nbin(1) - 1]
% disp("test_intersecting_bins_2d: For bin_idx 0, intersecting bins x: ");
% disp(bin_0_intersecting_x);
% unique_bin_0_intersecting_x = unique(bin_0_intersecting_x(bin_0_intersecting_x>=0 & bin_0_intersecting_x < grid_info.nbin(1)));
% expected_bin_0_intersecting_x = 0:(grid_info.nbin(1) - 1);
% disp("test_intersecting_bins_2d: For bin_idx 0, unique intersecting bins x: ");
% disp(unique_bin_0_intersecting_x);
% disp("test_intersecting_bins_2d: For bin_idx 0, expected intersecting bins x: ");
% disp(expected_bin_0_intersecting_x);
% assert(all(unique_bin_0_intersecting_x == expected_bin_0_intersecting_x));

% % Same for y
% unique_bin_0_intersecting_y = unique(bin_0_intersecting_y(bin_0_intersecting_y>=0 & bin_0_intersecting_y < grid_info.nbin(2)));
% expected_bin_0_intersecting_y = 0:(grid_info.nbin(2) - 1);
% assert(all(unique_bin_0_intersecting_y == expected_bin_0_intersecting_y));




%% test_0c
% Test that the returned bins for each query bin satisfy the proxy-circle
% Use grid_info and proxy_info from test_0b

grid_info_0c = grid_info;
proxy_info_0c = proxy_info;

N_bin_0c = grid_info.nbin(1) * grid_info.nbin(2);
all_bins = 0:(N_bin_0c - 1);
tol_dist = 1e-10;

for bin_idx = 0:(N_bin_0c - 1)
    [~, ~, binids] = intersecting_bins_2d(bin_idx, grid_info_0c, proxy_info_0c);
    c1 = bin_center(bin_idx, grid_info_0c);
    valid_returned = binids(binids >= 0);
    assert(length(valid_returned) == length(unique(valid_returned)), ...
        sprintf('bin_idx %d: returned bins have repeats', bin_idx));

    % First returned bins must intersect
    for k = 1:length(valid_returned)
        c2 = bin_center(valid_returned(k), grid_info_0c);
        d = norm(c1 - c2);
        assert(d <= 2 * proxy_info_0c.radius + tol_dist, ...
            sprintf('bin %d -> bin %d dist=%.6f > 2*r=%.6f', ...
                    bin_idx, valid_returned(k), d, 2*proxy_info_0c.radius));
    end

    %  non-returned valid bins must not intersect
    non_returned = setdiff(all_bins, valid_returned);
    for k = 1:length(non_returned)
        c2 = bin_center(non_returned(k), grid_info_0c);
        d = norm(c1 - c2);
        assert(d > 2 * proxy_info_0c.radius - tol_dist, ...
            sprintf('bin %d -> bin %d dist=%.6f <= 2*r=%.6f but not returned', ...
                    bin_idx, non_returned(k), d, 2*proxy_info_0c.radius));
    end
end


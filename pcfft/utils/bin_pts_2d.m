function [idstart,isort,r_sorted] = bin_pts_2d(r, dx, grid_info, nbin)
    % Sort points <r> into a number of bins.
    % nbin = number of gridpoints across that the bins are.
    if nargin < 4
        nbin = 1;
    end


    n_bins_per_dim = L * dx / nbin;

    % Find the ID of the bin in the X dim that each point occupies.
    id_x = round(r(1,:) / (nbin * dx)) + 1;
    id_y = round(r(2,:) / (nbin * dx)) + 1;

    bin_ids = id_x * n_bins_per_dim + id_y;
    [idx_bin_ids, sorted_bin_ids] = sort(bin_ids);

    % Sort the points
    r_sorted = r(:, idx_bin_ids);

    idbin = (id1-1)*nbin + id2;

    [idbinsort,isort] = sort(idbin);
    rsort = r(:,isort);

    % idstart = zeros(1,nbin*nbin+1);
    % ibin = 1; ibinold = 1; idstart(1) = 1;
    % for i = 1:size(r,2)
    %     if idbinsort(i)>ibin
    %         ibinold = ibin;
    %         ibin = idbinsort(i);
    %         idstart(ibinold+1:ibin) = i;
    %     end
    % end
    % idstart(ibin+1:end) = i+1;

end


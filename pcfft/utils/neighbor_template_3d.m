function [nbr_binids, nbr_gridpts, nbr_grididxes] = neighbor_template_3d(grid_info, proxy_info, bin_idx, template_pts, template_idxes)



    [~, ~, ~,nbr_binids] = intersecting_bins_3d(bin_idx, grid_info, proxy_info);
    % Compute the grid-index shift from center_bin to bin_idx.
    % cx, cy, cz are the 3D bin coordinates of center_bin
    cz = mod(grid_info.center_bin, grid_info.nbin(3));
    cy = mod(floor(grid_info.center_bin / grid_info.nbin(3)), grid_info.nbin(2));
    cx = floor(grid_info.center_bin / (grid_info.nbin(2) * grid_info.nbin(3)));


    % bz, by, bx are the 3D bin coordinates of bin_idx
    bz = mod(bin_idx, grid_info.nbin(3));
    by = mod(floor(bin_idx / grid_info.nbin(3)), grid_info.nbin(2));
    bx = floor(bin_idx / (grid_info.nbin(2) * grid_info.nbin(3)));

    delta_ix = (bx - cx) * grid_info.nbinpts;
    delta_iy = (by - cy) * grid_info.nbinpts;
    delta_iz = (bz - cz) * grid_info.nbinpts;

    shift = [delta_ix; delta_iy; delta_iz];
    shifted_idxes = (template_idxes.' + shift);
    % disp("neighbor_template_3d: size(template_idxes): " + int2str(size(template_idxes)) + ...
    %      ", size(shift): " + int2str(size(shift)) + ...
    %      ", size(shifted_idxes): " + int2str(size(shifted_idxes)));
    
    % Only keep in-bounds grid points.
    ngrid = grid_info.ngrid;
    in_bounds = shifted_idxes(1,:) >= 1 & shifted_idxes(1,:) <= ngrid(1) & ...
                shifted_idxes(2,:) >= 1 & shifted_idxes(2,:) <= ngrid(2) & ...
                shifted_idxes(3,:) >= 1 & shifted_idxes(3,:) <= ngrid(3);
    nbr_grididxes = (shifted_idxes(1,:) - 1) * ngrid(2) * ngrid(3) + (shifted_idxes(2,:) - 1) * ngrid(3) + shifted_idxes(3,:);
    % disp("neighbor_template_3d: size(nbr_grididxes): " + int2str(size(nbr_grididxes)) + ...
    %      ", size(in_bounds): " + int2str(size(in_bounds)));
    % Set out-of-bounds grid points to a dummy index.
    dummy_idx = ngrid(1) * ngrid(2) * ngrid(3) + 1;
    nbr_grididxes(~in_bounds) = dummy_idx;

    % Compute physical coordinates
    ctr = bin_center(bin_idx, grid_info);
    nbr_gridpts = template_pts + ctr;

    % disp("neighbor_template_3d: size(nbr_grididxes): " + int2str(size(nbr_grididxes)) + ...
        %  ", size(nbr_gridpts): " + int2str(size(nbr_gridpts)) + ...
        %  ", size(nbr_binids): " + int2str(size(nbr_binids)));

end
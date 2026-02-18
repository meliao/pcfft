classdef SortInfo
    %SORTINFO Information about sorting of source points into bins.
    % 
    % Sort points <src_info.r> into a number of bins.
    % Imagine a regular grid with bounds [xmin ymin xmax ymax] = Lbd
    % and grid spacing <dx>. There are [nx ny] = ngrid points 
    % in each dimension.
    % 
    % We want to sort the points into bins which are <nbinpts> 
    % regular gridpoints across. There are nbin(1) such bins in the x dimension
    % and nbin(2) in the y dimension.
    % Create indices for these bins by looping over x first and then y.
    %
    % Example of bin construction:
    % Suppose the points in r live on [-1, 1] x [-0.5, 0.5]
    % dx = 0.25
    % and if we set nbinpts = 3, we expect 
    % x bins [-1., -0.25], [-0.25, 0.5], [0.5, 1.]
    % y bins [-0.5, 0.25], [0.25, 0.5]
    %
    %
    %
    % Object properties:
    % r_srt: shape (2, n_src). Contains the source points sorted by bin
    %           index.
    % binid_srt: shape (n_src,). Contains the bin indexes of the sorted
    %           source points.
    % ptid_srt: shape (n_src,). Contains the indices of the sorting of 
    %           the source points. So r_srt = src_fin.r(:, ptid_srt).
    % id_start: shape (N_bins,). Contains the indices of the start of each
    %           bin in r_srt.
    % data_srt: struct of data associated with each point in r_srt. Each
    %           field must be of shape (l, n_srd) for some l
    %
    properties
        r_srt
        binid_srt
        ptid_srt
        id_start
        data_srt
    end

    methods
        function obj = SortInfo(src_info, dx, Lbd, nbin, nbinpts,der_fields)

            if nargin < 6 || isempty(der_fields), der_fields = {'r'}; end
            r = src_info.r;

            nbin = nbin(:);
            N_x_bins = nbin(1);
            N_y_bins = nbin(2);

            if size(r,1) == 2
                N_bins = N_x_bins * N_y_bins;


                % Find the ID of the bin in the X dim that each 
                % point occupies.
                % NOTE THAT id_x and id_y are zero-indexed!!
                id_x = floor((r(1,:) - Lbd(1)) / (nbinpts * dx));
                % Same for the Y dim.
                id_y = floor((r(2,:) - Lbd(2)) / (nbinpts * dx));

                bin_ids = id_x * N_y_bins + id_y;




            else
                % 3D code here.
                N_z_bins = nbin(3);
                N_bins = N_x_bins * N_y_bins * N_z_bins;

                % Find the ID of the bin in the X, Y, Z dims.
                % NOTE THAT id_x and id_y are zero-indexed!!
                id_x = floor((r(1,:) - Lbd(1)) / (nbinpts * dx));
                % Same for the Y dim.
                id_y = floor((r(2,:) - Lbd(2)) / (nbinpts * dx));
                % Same for the Z dim.
                id_z = floor((r(3,:) - Lbd(3)) / (nbinpts * dx));

                bin_ids = id_x * N_y_bins * N_z_bins + id_y * N_z_bins + id_z;
            end


            % Sort the bins
            [binid_srt, ptid_srt] = sort(bin_ids);

            % Sort the points
            r_srt = r(:, ptid_srt);

            % Form an array where id_start(i) gives us the index in 
            % r_sorted for the first point with bin idx i.
            % If the bin is empty, id_start(i) = id_start(i-1)
            % We want the slice (id_start(i+1) : id_start(i+2)-1) to give the
            % indices of points in bin i.
            id_start = ones(1,N_bins+1);

            % Loop through sorted_bin_ids and fill in id_start
            current_bin = 0;
            for i = 1:size(r,2)
                bin_i = binid_srt(i);
                if bin_i > current_bin
                    % Fill in all the bins we skipped
                    id_start(current_bin+2:bin_i+2) = i;
                    current_bin = bin_i;
                end
            end
            % Fill in the rest of the bins
            id_start(current_bin+2:N_bins+1) = size(r, 2) + 1;


            obj.r_srt = r_srt;
            obj.binid_srt = binid_srt;
            obj.ptid_srt = ptid_srt;
            obj.id_start = id_start;

            obj.data_srt = [];
            for field = der_fields
                obj.data_srt.(field{1}) = src_info.(field{1})(:,ptid_srt);
            end

        end
    end
end
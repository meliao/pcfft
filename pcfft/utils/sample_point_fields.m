function s = sample_point_fields(r, info, fields, idx)
    % Create synthetic points at R with fields copied from points IDX of INFO.
    %
    % Parameters
    % ----------
    % r : array [dim, n]
    %   positions of the synthetic points.
    % info : struct or object
    %   the real points (src_info or targ_info) to copy fields from.
    % fields : cell array of char
    %   names of the fields to copy. 'r' is ignored (it is always R).
    % idx : int array [1, n]
    %   which real points to copy each synthetic point's fields from.
    %
    % Returns
    % -------
    % s : struct
    %   s.r = r, and s.(f) is [ncomp, n] for every f in FIELDS.

    s = struct('r', r);
    if isempty(fields), return; end
    npts = size(info.r(:,:), 2);
    for k = 1:numel(fields)
        f = fields{k};
        if strcmp(f, 'r'), continue; end
        A = info.(f);
        A = A(:,:);
        if size(A, 2) == npts
            s.(f) = A(:, idx);
        elseif size(A, 1) == npts
            s.(f) = A(idx, :).';
        elseif isscalar(A)
            s.(f) = repmat(A, 1, numel(idx));
        else
            error('sample_point_fields: no dimension of field "%s" matches the number of points', f);
        end
    end
end

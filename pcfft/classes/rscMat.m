classdef rscMat
    % Custom sparse metrix format
    %
    % Attributes
    % ----------
    % 
    properties
        row_ptr
        col_ind
        vals
        nnz
        n
        m
    end

    methods
        function obj = rscMat(A)

            [i,j,vals] = find(A);
            [i,isort] = sort(i);
            j = j(isort);
            vals = vals(isort);
            obj.nnz = length(vals);
            obj.row_ptr = zeros(1,size(A,1)+1);
            % idiff = diff(i);
            % ijumps = find(idiff);
            obj.row_ptr(1) = 1;
            ctr = 1;
            for l = 1:size(A,1)
                % if i(ctr+1) > l
                %     obj.row_ptr(l+1) = obj.row_ptr(l);
                % else
                ni = 0;
                for k = ctr:obj.nnz
                    if i(k) > l, break; end
                    ni = ni + 1;
                end
                obj.row_ptr(l+1) = obj.row_ptr(l) + ni;
                ctr = ctr + ni;
            end
            % ctr = 1;
            % for i = 1:length(ijumps)
            %     nis = ijumps(i+1)-ijumps(i);
            %     obj.row_ptr(i+1) = obj.row_ptr(i+1);
            %     ctr = nis + ctr;
            % end
            obj.col_ind = j;
            
            obj.vals = vals;
            obj.m = size(A,1);
            obj.n = size(A,2);
        end
        function a = mtimes(obj,b)
            % a = zeros(obj.m,1);
            % for l = 1:obj.m
            %     ivals = obj.row_ptr(l):obj.row_ptr(l+1)-1;
            %     a(l) = obj.vals(ivals)*b(obj.col_ind(ivals));
            % end

            % c = obj.vals.*b(obj.col_ind).';
            % a = zeros(obj.m,1);
            % for l = 1:obj.m
            %     ivals = obj.row_ptr(l):obj.row_ptr(l+1)-1;
            %     a(l) = sum(c(ivals));
            % end

            c = obj.vals.*b(obj.col_ind);
            a = zeros(obj.m,1);
            for l = 1:obj.m
                a(l) = sum(c(obj.row_ptr(l):obj.row_ptr(l+1)-1));
            end
        end
    end
end

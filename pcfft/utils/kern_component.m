function val = kern_component(kern, src, targ,i,j,opdim)
    % This routine evaluates the ijth component of a matrix valued kernel.
    %      i.e. K_{i,j}(src, targ)
    % It assumes that the values are kern values i.e.
    % the submatrix val(1:opdim(1),1:opdim(2)) gives to the interaction
    % between the first source-target pairs
    %
    % Parameters
    % ----------
    % kern : kernel
    %   The matrix valued kernel.
    % src : struct
    %   Specifies the source points.
    % targ : struct
    %   Specifies the source points.
    % i : integer
    %   the target component
    % j : integer
    %   the source component to take
    % opdim : integer array [1,2]
    %   dimension of the kernel matrix evaluated at each source-target pair
    %
    % Returns
    % -------
    % val : matrix [ntarg, nsrc]
    %   
    if ~isa(kern,'function_handle')
        try
            kern = kern.eval;
        catch
            error('kern is not a function and does not have an eval property')
        end
    end
    val = kern(src,targ);

    val = val(i:opdim(1):end, j:opdim(2):end);

end
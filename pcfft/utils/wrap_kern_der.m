function val = wrap_kern_der(kern, src, targ, nder, h)
    % This routine evaluates the kernel and its first nder derivatives. The
    % derivatives are computed
    %
    % Parameters
    % ----------
    % kern : kernel
    %   The matrix valued kernel.
    % src : struct
    %   Specifies the source points.
    % targ : struct
    %   Specifies the source points.
    % nder : int
    %   Number of derivatives to return
    % h : float (optional)
    %   finite difference step size. Defaults to 1e-2.
    %
    % Returns
    % -------
    % val : matrix [opdims(1)*ntarg, opdims(2)*nsrc]
    %   
    if nargin < 5
        h = 1e-2;
    end
    if ~isa(kern,'function_handle')
        try
            kern = kern.eval;
        catch
            error('kern is not a function and does not have an eval property')
        end
    end
    val0 = kern(src,targ);

    if nder > 0
        ntarg = targ.r(:,:)./vecnorm(targ.r(:,:));
        
        targp = targ;
        targp.r = targ.r + h*ntarg;
        
        targm = targ;
        targm.r = targ.r - h*ntarg;
        
        valp = kern(src,targp);
        valm = kern(src,targm);
    end

    if nder == 0
        val = val0;
    elseif nder == 1
        val = zeros(2*size(val0,1), size(val0,2));
        val(1:2:end,:) = val0;
        val(2:2:end,:) = (valp-valm) / 2 /h;
    elseif nder == 2
        val = zeros(3*size(val0,1), size(val0,2));
        val(1:3:end,:) = val0;
        val(2:3:end,:) = (valp-valm) / 2 /h;
        val(3:3:end,:) = (valp-2*val0+valm) / h^2;
    else
        error('wrap_kern_der: Too many proxy derivatives')
    end



end
function [idstart,isort,rsort] = bin_pts(r,rbin,dx,nbin,L)

    id1 = round((r(1,:)+L)/(rbin*dx))+1;
    id2 = round((r(2,:)+L)/(rbin*dx))+1;


    idbin = (id1-1)*nbin + id2;

    [idbinsort,isort] = sort(idbin);
    rsort = r(:,isort);

    idstart = zeros(1,nbin*nbin+1);
    ibin = 1; ibinold = 1; idstart(1) = 1;
    for i = 1:size(r,2)
        if idbinsort(i)>ibin
            ibinold = ibin;
            ibin = idbinsort(i);
            idstart(ibinold+1:ibin) = i;
        end
    end
    idstart(ibin+1:end) = i+1;

end


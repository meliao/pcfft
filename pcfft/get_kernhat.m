function kern_hat = get_kernhat(kern,rgrid, ngrid, Lbd, dx)
% evaluate Fourier transform of kernel on rgrid

rgrid0 = Lbd(:,1) + dx*ngrid(:);
kernvals = kern(rgrid0, rgrid);

kernvals = reshape(kernvals, flip(2*ngrid(:)'+1));

kernvalshift = ifftshift(kernvals);
kern_hat = fftn(kernvalshift);

end
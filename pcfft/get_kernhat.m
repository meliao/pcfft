function kern_hat = get_kernhat(kern,rgrid, ngrid)
% evaluate Fourier transform of kernel on grid

kernvals = kern(rgrid(:,1), rgrid);
kernvals = reshape(kernvals, ngrid,ngrid);

kernvalshift = fftshift(kernvals);
kern_hat = fft2(kernvalshift);

end
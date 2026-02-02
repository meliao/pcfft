function kern_hat = get_kernhat(kern_0, ngrid, Lbd, dx, offset)
% evaluate Fourier transform of free-space kernel on rgrid

rgrid0 = Lbd(:,1) + dx*ngrid(:) - offset;

% get padded grid
if length(ngrid) == 2
    % Create a regular grid with spacing dx starting at the xmin, ymin point
    % specified by Lbd. 
    xx = Lbd(1, 1) - offset + (0: 2*ngrid(1) - 1) * dx;
    yy = Lbd(2, 1) - offset + (0: 2*ngrid(2) - 1) * dx;
    [X, Y] = meshgrid(xx, yy);
    rgrid = [X(:).'; Y(:).'];
elseif length(ngrid) == 3
    xx = Lbd(1, 1) - offset + (0: 2*ngrid(1) - 1) * dx;
    yy = Lbd(2, 1) - offset + (0: 2*ngrid(2) - 1) * dx;
    zz = Lbd(3, 1) - offset + (0: 2*ngrid(3) - 1) * dx;
    [X, Y, Z] = meshgrid(xx, yy, zz);
    X = permute(X,[3,1,2]);
    Y = permute(Y,[3,1,2]);
    Z = permute(Z,[3,1,2]);
    rgrid = [X(:).'; Y(:).'; Z(:).'];
end
% evaluate kernel
kernvals = kern_0(struct('r',rgrid0), struct('r',rgrid));

% Fourier transform
kernvals = reshape(kernvals, 2*flip(ngrid(:)'));
kernvals = ifftshift(kernvals);
kern_hat = fftn(kernvals);

end
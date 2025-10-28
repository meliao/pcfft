function [pxy,str2pxy,pxygrid] = setup_proxy(dxgrid,npxy,rpxy,gkern,tol,zk)
if nargin < 5
    tol = 1e-10;
end
if nargin < 6
    zk = [];
end
if isempty(zk)
    pxy = get_proxy_ann(rpxy,20,51);
    % pxy = get_proxy_ann(rpxy,20,31);
    % pxy = get_proxy_ann(rpxy,40,51);
else
    % pxy = get_proxy_ann(rpxy,10,11,max(0.5*2*pi/zk/rpxy,1.5));
    % pxy = get_proxy_ann(rpxy,10,21,max(0.5*2*pi/zk/rpxy,1.5));
    % pxy = get_proxy_ann(rpxy,40,91,max(2*pi/zk/rpxy,1.5));
    pxy = get_proxy_ann(rpxy,5,91,max(2*pi/zk/rpxy,1.5));
end

[pxygrid_x,pxygrid_y] = meshgrid(linspace(-(npxy-1)/2*dxgrid,(npxy-1)/2*dxgrid,npxy));
pxygrid_z = 0*pxygrid_x;
pxygrid = [pxygrid_x(:).';pxygrid_y(:).'; pxygrid_z(:).'];

grid2pxy = gkern.eval(struct('r', pxygrid), struct('r',pxy));

% [U,sigma,V] = svd(grid2pxy,"econ");
% singval = diag(sigma); 
% singval(singval<tol) = Inf; % hacky way to zero out small singular values
% sigma(1:size(sigma,1)+1:end) = 1./singval;
% pxy2str = V*sigma*U';

% [sk,rd,T] = id(grid2pxy.',tol);
sk = 1:size(pxy,2);
% norm(grid2pxy(rd,:) - T.'*grid2pxy(sk,:))
% [Q,R] = qr(grid2pxy(sk,:));
% 
% % pxy2pxyskel = eye(size(pxy,2)); pxy2pxyskel = pxy2pxyskel(sk,:);
% % pxy2str = R\(Q'*pxy2pxyskel);
% 

% pxy2str = R\(Q');
pxy = pxy(:,sk);
str2pxy = grid2pxy(sk,:);

% %%% Least squares test
% src = [dxgrid/3;dxgrid/2;0];
% src2pxy = gkern.eval(struct('r', src), struct('r',pxy));
% 
% targ = rpxy*2 * [cos(0.1); sin(0.1); 0];
% 
% src2targ = gkern.eval(struct('r', src), struct('r',targ));
% grid2targ = gkern.eval(struct('r', pxygrid), struct('r',targ));
% norm(src2targ - grid2targ*pxy2str*src2pxy)
% norm(src2targ - grid2targ*(grid2pxy\src2pxy))
end

function pxy = get_proxy_ann(rmin,nrad,ntheta,rmax)
    if nargin < 2
        nrad = 6;
    end
    if nargin < 3
        ntheta = 31;
        ntheta = 51;
    end
    if nargin < 4
        rmax = 1.5;
        % rmax = 5;
    end

    thetas = linspace(0,2*pi,ntheta+1); thetas = thetas(1:end-1);
    % rads = linspace(rmin, rmin*rmax,nrad).';
    rads = 1./linspace(1/(rmin*rmax), 1/rmin,nrad).';
    
    pxy_x = rads .* cos(thetas);
    pxy_y = rads .* sin(thetas);
    pxy = [pxy_x(:).'; pxy_y(:).';0*pxy_x(:).'];
end

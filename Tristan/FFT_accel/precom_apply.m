function [Aspread,Asubtract,kern_hat,rsort,isort,times] = precom_apply(gkern,zk,S,L,dx,ngrid,rgrid,npxy,rpxy,rbin,iftimes,tol)
if nargin < nargin('precom_apply')-1
    iftimes = 0;
end

if nargin < nargin('precom_apply')
    tol = 1e-12;
    tol = 1e-13;
end

tic;
[pxy,str2pxy,pxygrid] = setup_proxy(dx,npxy,rpxy,gkern,tol,zk);
tproxy = toc;


nbin = floor(ngrid/rbin);

tic;
[idstart,isort,rsort] = bin_pts(S.r,rbin,dx,nbin,L);


truebin = floor(1*nbin/6):ceil(5*nbin/6);
bin_ctrs = zeros(3,nbin*nbin);
for i = 1:nbin
    xcenter = rbin * (i-1) * dx - L;
    for j = 1:nbin
        ycenter =  rbin * (j-1) * dx - L;
        rcenter = [xcenter;ycenter;0];
        bin_ctrs(:,(i-1)*nbin + j)= rcenter;
    end
end
tbin = toc;

tic;
r_loc = 0*rsort;
for ibin = 1:nbin*nbin
    if idstart(ibin) == idstart(ibin+1), continue, end

    id_loc = idstart(ibin):(idstart(ibin+1)-1);
    r_loc(:,id_loc) = rsort(:,id_loc) - bin_ctrs(:,ibin);
end

pts2pxy = gkern(struct('r',r_loc),struct('r',pxy));

% pts2grid = pxy2str*pts2pxy;
pts2grid = str2pxy\pts2pxy;

Aspread = sparse(ngrid*ngrid,S.npts);
for id_bin_x = truebin
    for id_bin_y = truebin
    ibin = (id_bin_x-1)*nbin + id_bin_y;
    if idstart(ibin) == idstart(ibin+1), continue, end

    id_loc = idstart(ibin):(idstart(ibin+1)-1);

    id_grid_x = (id_bin_x-1)*rbin + (-((npxy-1)/2):((npxy-1)/2))+1;
    id_grid_y = (id_bin_y-1)*rbin + (-((npxy-1)/2):((npxy-1)/2))+1;

    id_grid_loc = (id_grid_x(:).'-1)*ngrid + id_grid_y(:);

    Aspread(id_grid_loc(:),id_loc) = Aspread(id_grid_loc(:),id_loc) + pts2grid(:,id_loc);
    end
end

tpre = toc;


tic;
kernvals = gkern(struct('r',[0;0;0]), struct('r',rgrid));
kernvals = reshape(kernvals, ngrid,ngrid);

kernvalshift = fftshift(kernvals);
kern_hat = fft2(kernvalshift);
tkern = toc;

nsub_rad = ceil(rpxy/dx);

nsub_bin = ceil((nsub_rad+(npxy-1)/2) / rbin);

nloc_grid = (nsub_bin*rbin + (npxy-1)/2);
[rgrid_sub_x,rgrid_sub_y] = meshgrid((-nloc_grid:nloc_grid)*dx);
rgrid_sub = [rgrid_sub_x(:).';rgrid_sub_y(:).';0*rgrid_sub_x(:).'];

[rbin_sub_x,rbin_sub_y] = meshgrid((-nsub_bin:nsub_bin)*rbin*dx);
rbin_sub = [rbin_sub_x(:).';rbin_sub_y(:).';0*rbin_sub_x(:).'];

% figure(2);clf
% my_scatter(rgrid_sub)
% hold on
% my_scatter(pxygrid,'.','linewidth',2)
% my_scatter(pxy,'*','linewidth',2)
% my_scatter(rbin_sub,'x','linewidth',2)
% 
% my_scatter(pxy+[nsub_bin*rbin*dx;0;0],'*','linewidth',2)
% my_scatter(pxygrid+[nsub_bin*rbin*dx;0;0],'o','linewidth',2)
% hold off

t1 = tic;
loc_kern = gkern(struct('r',pxygrid),struct('r',rgrid_sub));

idnear_startbin = ones(1,nbin*nbin); idnear_startbin(1) = 1;
ilasttouch = 0;
rnear = [];
for id_bin_x = truebin
    for id_bin_y = truebin
        ibin = (id_bin_x-1)*nbin + id_bin_y;
        
        if idstart(ibin) == idstart(ibin+1)
            if ibin == ilasttouch+1
                % updating the rest of idstart everytime is slow
            idnear_startbin(ibin+1:end) = idnear_startbin(ibin);
            end
            continue
        end

        id_loc = idstart(ibin):(idstart(ibin+1)-1);
    
        rsrc = rsort(:,id_loc);

        ibins = (id_bin_x+(-nsub_bin:nsub_bin)-1)*nbin + id_bin_y+(-nsub_bin:nsub_bin).';
        ibins = ibins(:);
        idstart_loc = zeros(1,length(ibins)+1);
        idstart_loc(1) = 1;
        id_targs = [];
        for i = 1:length(ibins)
            idstart_loc(i+1) = idstart_loc(i) + idstart(ibins(i)+1)-idstart(ibins(i));
            id_targs = [id_targs,idstart(ibins(i)):(idstart(ibins(i)+1)-1)];
        end
        rnear_loc = rsort(:,id_targs) - reshape(rsrc,3,1,[]);
        rnear = [rnear, rnear_loc(:,:)];
        idnear_startbin(ibin+1) = idnear_startbin(ibin) + size(rnear_loc(:,:),2);
        
        ilasttouch =ibin;
    end
end

sig2u = gkern(struct('r',[0;0;0]),struct('r',rnear));

iid = zeros(1,size(rnear,2));
jid = zeros(1,size(rnear,2));
vals = zeros(1,size(rnear,2));


pts2subpts = loc_kern*pts2grid;
for id_bin_x = truebin
    for id_bin_y = truebin
        ibin = (id_bin_x-1)*nbin + id_bin_y;
        if idstart(ibin) == idstart(ibin+1), continue, end

        id_loc = idstart(ibin):(idstart(ibin+1)-1);
        % sig2loc = loc_kern*pts2grid(:,id_loc);
        sig2loc = pts2subpts(:,id_loc);
        % a = pts2grid(id_loc,:).';
    
        % rsrc = rsort(:,id_loc);

        ibins = (id_bin_x+(-nsub_bin:nsub_bin)-1)*nbin + id_bin_y+(-nsub_bin:nsub_bin).';
        ibins = ibins(:);
        idstart_loc = zeros(1,length(ibins)+1);
        idstart_loc(1) = 1;
        id_targs = [];
        for i = 1:length(ibins)
            idstart_loc(i+1) = idstart_loc(i) + idstart(ibins(i)+1)-idstart(ibins(i));
            id_targs = [id_targs,idstart(ibins(i)):(idstart(ibins(i)+1)-1)];
        end

        id_sparse = idnear_startbin(ibin):(idnear_startbin(ibin+1)-1);
        id_sparse = reshape(id_sparse, [], length(id_loc));
        sig2ubin = sig2u(id_sparse);

        Aloc2 = zeros(size(id_sparse));
        
        for i = -nsub_bin:nsub_bin
            id_grid_x = i*rbin + (-((npxy-1)/2):((npxy-1)/2))+1+nloc_grid;
            
            for j = -nsub_bin:nsub_bin

            ibin = (id_bin_x+i-1)*nbin + id_bin_y+j;
            ibinloc = (i+nsub_bin+1-1)*(2*nsub_bin+1) + j+nsub_bin+1;
            ipt_loc = idstart_loc(ibinloc):(idstart_loc(ibinloc+1)-1);
            id_loc_targ = idstart(ibin):(idstart(ibin+1)-1);

            
            id_grid_y = j*rbin + (-((npxy-1)/2):((npxy-1)/2))+1+nloc_grid;
            id_grid_loc = (id_grid_x(:).'-1)*(2*nloc_grid+1) + id_grid_y(:);

            sig2targ = pts2grid(:,id_loc_targ).' * sig2loc(id_grid_loc(:),:);
                        
            sig2uloc = sig2ubin(ipt_loc,:);

            Aloc2(ipt_loc,:) = sig2uloc - sig2targ;
            end
        end
        
        is = repmat(id_targs(:), 1, length(id_loc));
        js = repmat(id_loc(:).', length(id_targs), 1);
        iid(id_sparse) = is(:).';
        jid(id_sparse) = js(:).';
        vals(id_sparse) = Aloc2(:).';
    end
end
Asubtract = sparse(iid, jid, vals);

% toc;
tsub_pre = toc(t1);


times = [tsub_pre,tpre,tproxy,tbin,tkern];

if iftimes
    % %%
    % fprintf('add and subtract op =  %.2e s\n',tsub_pre)
    % fprintf('spread op = %.2e s\n',tpre)
    % fprintf('proxy = %.2e s\n',tproxy)
    % fprintf('binning = %.2e s\n',tbin)
    % fprintf('fft kern =  %.2e s\n',tkern)
    
    % tproxy
    % tbin
    tkern
    tpre
    tsub_pre
    % %%
end

end
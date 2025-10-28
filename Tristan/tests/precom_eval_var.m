function [spread_info,times] = precom_eval_var(gkern_0,gkern,zk,srcinfo,targinfo,L,dx,ngrid,rgrid,npxy,rpxy,rbin,ib2v,iftimes,tol)
% function [Aspread_s,Aspread_t,Asubtract,kern_hat,src_sort,isort,targ_sort,isort_t,times] = precom_eval_var(gkern_0,zk,srcinfo,gkern_s,targinfo,gkern_t,gkern_st,L,dx,ngrid,rgrid,npxy,rpxy,rbin,ib2v,iftimes,tol)
if nargin < nargin('precom_eval_var')-1
    ib2v = 0;
end
if nargin < nargin('precom_eval_var')-1
    iftimes = 0;
end

if nargin < nargin('precom_eval_var')
    tol = 1e-12;
    tol = 1e-8;
end

gkern = kernel(gkern);
gkern_0 = kernel(gkern_0);

if ib2v
    gkern_s = gkern;
    gkern_t = gkern_0;
    gkern_st = gkern;
else
    gkern_s = gkern_0;
    gkern_r = kernel(@(s,t) gkern.eval(t,s).');
    gkern_t = gkern_r;
    gkern_st = gkern_r;
end

  warnStruct = warning('off','MATLAB:rankDeficientMatrix');


tic;
[pxy,str2pxy,pxygrid] = setup_proxy(dx,npxy,rpxy,gkern_0,tol,zk);
tproxy = toc;


nbin = floor(ngrid/rbin);
%tic;
[Aspread_s,pts2grid,idstart,isort,src_sort,truebin,tbin,tpre] = getspread(gkern_s,srcinfo,rbin,dx,nbin,L,pxy,str2pxy,ngrid,npxy);
%toc;

%%
[Aspread_t,pts2grid_t,idstart_t,isort_t,targ_sort,truebin,tbin_t,tpre_t] = getspread(gkern_t,targinfo,rbin,dx,nbin,L,pxy,str2pxy,ngrid,npxy);

tic;
kernvals = gkern_0.eval(struct('r',[0;0;0]), struct('r',rgrid));
kernvals = reshape(kernvals, ngrid,ngrid);

kernvalshift = fftshift(kernvals);
kern_hat = fft2(kernvalshift);
tkern = toc;


if ib2v
    [Asubtract,tsub_pre] = get_add_sub(gkern_0,gkern_st,pxygrid,rpxy,dx,rbin,npxy,nbin,idstart,src_sort,idstart_t,targ_sort,pts2grid,pts2grid_t,truebin,ib2v);
else
    [Asubtract,tsub_pre] = get_add_sub(gkern_0,gkern_st,pxygrid,rpxy,dx,rbin,npxy,nbin,idstart_t,targ_sort,idstart,src_sort,pts2grid_t,pts2grid,truebin,ib2v);
    Asubtract = Asubtract.';
end

  warning(warnStruct);

spread_info = [];
spread_info.Aspread_s =  Aspread_s;
spread_info.Aspread_t =  Aspread_t;
spread_info.Asubtract =  Asubtract;
spread_info.kern_hat =  kern_hat;
spread_info.src_sort =  src_sort;
spread_info.isort =  isort;
spread_info.targ_sort =  targ_sort;
spread_info.isort_t =  isort_t;
spread_info.ngrid = ngrid;

times = [tsub_pre,tpre,tproxy,tbin,tkern];

if iftimes
    % %%
    fprintf('fft kern =%.2e s\n',tkern)
    fprintf('spread op = %.2e s\n',tpre)
    fprintf('add and subtract op = %.2e s\n',tsub_pre)
    % fprintf('proxy = %.2e s\n',tproxy)
    % fprintf('binning = %.2e s\n',tbin)
    
    
    % tproxy
    % tbin
    % tkern
    % tpre
    % tsub_pre
    % %%
end

end


function [Aspread,pts2grid,idstart,isort,pt_sort,truebin,tbin,tpre] = getspread(gkern_p,ptinfo,rbin,dx,nbin,L,pxy,str2pxy,ngrid,npxy)

tic;
[idstart,isort,rsort] = bin_pts(ptinfo.r(:,:),rbin,dx,nbin,L);
pt_sort = [];
pt_sort.r = ptinfo.r(:,isort);
if isfield(ptinfo,'n') || isprop(ptinfo,'n')
    pt_sort.n = ptinfo.n(:,isort);
end
if isfield(ptinfo,'d') || isprop(ptinfo,'d')
    pt_sort.d = ptinfo.d(:,isort);
end
if isfield(ptinfo,'d2') || isprop(ptinfo,'d2')
    pt_sort.d2 = ptinfo.d2(:,isort);
end

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
pt_loc = [];
pt_loc.r = pt_sort.r;
if isfield(pt_sort,'n') 
    pt_loc.n = pt_sort.n;
end
if isfield(pt_sort,'d') 
    pt_loc.d = pt_sort.d;
end
if isfield(pt_sort,'d2') 
    pt_loc.d2 = pt_sort.d2;
end
for ibin = 1:nbin*nbin
    if idstart(ibin) == idstart(ibin+1), continue, end

    id_loc = idstart(ibin):(idstart(ibin+1)-1);
    pt_loc.r(:,id_loc) = pt_sort.r(:,id_loc) - bin_ctrs(:,ibin);
end

[mm, nn] = size(str2pxy);
npt = size(pt_loc.r(:,:),2);
pts2grid = zeros(nn, npt);
pxy_struct = [];
pxy_struct.r = pxy;
nbsize = 100;
nbat = ceil(npt/nbsize);
for i = 1:nbat
    istart = (i-1)*nbsize+1;
    iend = min(i*nbsize,npt);
    iind = istart:iend;
    % fprintf('istart=%d   iend=%d\n',istart,iend);
    pstruct = [];
    pstruct.r = pt_loc.r(:,iind);

    if isfield(pt_loc,'n') 
        pstruct.n = pt_loc.n(:,iind);
    end
    if isfield(pt_loc,'d') 
        pstruct.d = pt_loc.d(:,iind);
    end
    if isfield(pt_loc,'d2') 
        pstruct.d2 = pt_loc.d2(:,iind);
    end

    ptmp = gkern_p.eval(pstruct, pxy_struct);
    pts2grid(:,iind) = str2pxy\ptmp;
end
% pts2pxy = gkern_p.eval(pt_loc,struct('r',pxy));

% pts2grid = pxy2str*pts2pxy;
% pts2grid = str2pxy\pts2pxy;

Aspread = spalloc(ngrid*ngrid,size(ptinfo.r(:,:),2),npxy^2*size(ptinfo.r(:,:),2));
%I = [];
%J = [];
%v = [];
% Aspread = sparse(ngrid*ngrid,size(ptinfo.r(:,:),2));
Aspread = spalloc(ngrid*ngrid,size(ptinfo.r(:,:),2),npxy^2*size(ptinfo.r(:,:),2));
for id_bin_x = truebin
    for id_bin_y = truebin
    ibin = (id_bin_x-1)*nbin + id_bin_y;
    if idstart(ibin) == idstart(ibin+1), continue, end

    id_loc = idstart(ibin):(idstart(ibin+1)-1);

    id_grid_x = (id_bin_x-1)*rbin + (-((npxy-1)/2):((npxy-1)/2))+1;
    id_grid_y = (id_bin_y-1)*rbin + (-((npxy-1)/2):((npxy-1)/2))+1;

    id_grid_loc = (id_grid_x(:).'-1)*ngrid + id_grid_y(:);
%    [jj, ii] = meshgrid(id_loc(:),id_grid_loc(:));
%    vals = pts2grid(:,id_loc);

%    I = [I; ii(:)];
%    J = [J; jj(:)];
    
%    v = [v; vals(:)];
    
    Aspread(id_grid_loc(:),id_loc) = pts2grid(:,id_loc);
    end
end
% Aspread = sparse(I,J,v,ngrid*ngrid,size(ptinfo.r(:,:),2));
tpre = toc;

end


function [Asubtract,tsub_pre] = get_add_sub(gkern_0,gkern_st,pxygrid,rpxy,dx,rbin,npxy,nbin,idstart,src_sort,idstart_t,targ_sort,pts2grid,pts2grid_t,truebin,ib2v)
ib2v = 1;
nsub_rad = ceil(rpxy/dx);

nsub_bin = ceil((nsub_rad+(npxy-1)/2) / rbin);

nloc_grid = (nsub_bin*rbin + (npxy-1)/2);
[rgrid_sub_x,rgrid_sub_y] = meshgrid((-nloc_grid:nloc_grid)*dx);
rgrid_sub = [rgrid_sub_x(:).';rgrid_sub_y(:).';0*rgrid_sub_x(:).'];

[rbin_sub_x,rbin_sub_y] = meshgrid((-nsub_bin:nsub_bin)*rbin*dx);
rbin_sub = [rbin_sub_x(:).';rbin_sub_y(:).';0*rbin_sub_x(:).'];

% figure(1);clf
% my_scatter(rgrid_sub)
% hold onsr
% my_scatter(pxygrid,'.','linewidth',2)
% my_scatter(pxy,'*','linewidth',2)
% my_scatter(rbin_sub,'x','linewidth',2)
% 
% my_scatter(pxy+[nsub_bin*rbin*dx;0;0],'*','linewidth',2)
% my_scatter(pxygrid+[nsub_bin*rbin*dx;0;0],'o','linewidth',2)
% hold off

t1 = tic;
loc_kern = gkern_0.eval(struct('r',pxygrid),struct('r',rgrid_sub));

idnear_startbin = ones(1,nbin*nbin); idnear_startbin(1) = 1;
% id_startbin_t = ones(1,nbin*nbin); id_startbin_t(1) = 1;
ilasttouch = 0;
rnear = [];
nnear = [];
dnear = [];
d2near = [];

for id_bin_x = truebin
    for id_bin_y = truebin
        ibin = (id_bin_x-1)*nbin + id_bin_y;
        
        if (idstart(ibin) == idstart(ibin+1)) %|| (idstart_t(ibin) == idstart_t(ibin+1))
            if ibin == ilasttouch+1
                % updating the rest of idstart everytime is slow
            idnear_startbin(ibin+1:end) = idnear_startbin(ibin);
            % id_startbin_t(ibin+1:end) = id_startbin_t(ibin);
            end
            % idnear_startbin(ibin+1) = idnear_startbin(ilasttouch+1);
            continue
        end

        id_loc = idstart(ibin):(idstart(ibin+1)-1);
    
        rsrc_loc = src_sort.r(:,id_loc);

        ibins = (id_bin_x+(-nsub_bin:nsub_bin)-1)*nbin + id_bin_y+(-nsub_bin:nsub_bin).';
        ibins = ibins(:);
        idstart_loc_s = zeros(1,length(ibins)+1); idstart_loc_s(1) = 1;
        idstart_loc_t = zeros(1,length(ibins)+1); idstart_loc_t(1) = 1;
        id_targs = [];
        for i = 1:length(ibins)
            idstart_loc_s(i+1) = idstart_loc_s(i) + idstart(ibins(i)+1)-idstart(ibins(i));
            idstart_loc_t(i+1) = idstart_loc_t(i) + idstart_t(ibins(i)+1)-idstart_t(ibins(i));
            id_targs = [id_targs,idstart_t(ibins(i)):(idstart_t(ibins(i)+1)-1)];
        end
        rnear_loc = targ_sort.r(:,id_targs) - reshape(rsrc_loc,3,1,[]);
        rnear = [rnear, rnear_loc(:,:)];
        if ~ib2v 
            if isfield(targ_sort,'n')
                nnear_loc = repmat(targ_sort.n(:,id_targs),1,1,length(id_loc));
                nnear = [nnear, nnear_loc(:,:)];
            end
            if isfield(targ_sort,'d')
                dnear_loc = repmat(targ_sort.d(:,id_targs),1,1,length(id_loc));
                dnear = [dnear, dnear_loc(:,:)];
            end
            if isfield(targ_sort,'d2')
                d2near_loc = repmat(targ_sort.d2(:,id_targs),1,1,length(id_loc));
                d2near = [d2near, d2near_loc(:,:)];
            end
        else
            if isfield(src_sort,'n')
                nnear_loc = reshape(src_sort.n(:,id_loc),[],1,length(id_loc));
                nnear_loc = repmat(nnear_loc,1,length(id_targs),1);
                nnear = [nnear, nnear_loc(:,:)];
            end
            if isfield(src_sort,'d')
                dnear_loc = reshape(src_sort.d(:,id_loc),[],1,length(id_loc));
                dnear_loc = repmat(dnear_loc,1,length(id_targs),1);
                dnear = [dnear, dnear_loc(:,:)];
            end
            if isfield(src_sort,'d2')
                d2near_loc = reshape(src_sort.d2(:,id_loc),[],1,length(id_loc));
                d2near_loc = repmat(d2near_loc,1,length(id_targs),1);
                d2near = [d2near, d2near_loc(:,:)];
            end
        end
        idnear_startbin(ibin+1) = idnear_startbin(ibin) + size(rnear_loc(:,:),2);
        
        ilasttouch =ibin;
    end
end

nearinfo = [];
nearinfo.r = rnear;
% nearinfo.n = nnear; % added here by me
if ~isempty(nnear)
    nearinfo.n = nnear;
end
if ~isempty(dnear)
    nearinfo.d = dnear;
end
if ~isempty(d2near)
    nearinfo.d2 = d2near;
end

npt = size(nearinfo.r(:,:),2);
sig2u = zeros(npt,1);
nbsize = 10000;
nbat = ceil(npt/nbsize);

if ib2v && isfield(src_sort,'n')
    nearinfo.r = -nearinfo.r;  
    for i = 1:nbat
        istart = (i-1)*nbsize+1;
        iend = min(i*nbsize,npt);
        iind = istart:iend;
        % fprintf('istart=%d   iend=%d\n',istart,iend);
        pstruct = [];
        pstruct.r = nearinfo.r(:,iind);
        if ~isempty(nnear)
            pstruct.n = nearinfo.n(:,iind);
        end
        if ~isempty(dnear)
            pstruct.d = nearinfo.d(:,iind);
        end
        if ~isempty(d2near)
            pstruct.d2 = nearinfo.d2(:,iind);
        end
        sig2u(iind) = gkern_st.eval(pstruct, struct('r',[0;0;0])).';
    end
else
    for i = 1:nbat
        istart = (i-1)*nbsize+1;
        iend = min(i*nbsize,npt);
        iind = istart:iend;
        % fprintf('istart=%d   iend=%d\n',istart,iend);
        pstruct = [];
        pstruct.r = nearinfo.r(:,iind);
        if ~isempty(nnear)
            pstruct.n = nearinfo.n(:,iind);
        end
        if ~isempty(dnear)
            pstruct.d = nearinfo.d(:,iind);
        end
        if ~isempty(d2near)
            pstruct.d2 = nearinfo.d2(:,iind);
        end
        sig2u(iind) = gkern_st.eval(struct('r',[0;0;0]),pstruct);
    end
end

iid = zeros(1,size(rnear,2));
jid = zeros(1,size(rnear,2));
vals = zeros(1,size(rnear,2));


pts2subpts = loc_kern*pts2grid;
for id_bin_x = truebin
    for id_bin_y = truebin
        ibin = (id_bin_x-1)*nbin + id_bin_y;
        if (idstart(ibin) == idstart(ibin+1))% || (idstart_t(ibin) == idstart_t(ibin+1))
            continue
        end

        id_loc = idstart(ibin):(idstart(ibin+1)-1);
        % sig2loc = loc_kern*pts2grid(:,id_loc);
        sig2loc = pts2subpts(:,id_loc);
        % a = pts2grid(id_loc,:).';
    
        % rsrc = rsort(:,id_loc);

        ibins = (id_bin_x+(-nsub_bin:nsub_bin)-1)*nbin + id_bin_y+(-nsub_bin:nsub_bin).';
        ibins = ibins(:);
        idstart_loc_s = zeros(1,length(ibins)+1); idstart_loc_s(1) = 1;
        idstart_loc_t = zeros(1,length(ibins)+1); idstart_loc_t(1) = 1;
        for i = 1:length(ibins)
            idstart_loc_s(i+1) = idstart_loc_s(i) + idstart(ibins(i)+1)-idstart(ibins(i));
            idstart_loc_t(i+1) = idstart_loc_t(i) + idstart_t(ibins(i)+1)-idstart_t(ibins(i));
        end

        id_targs = zeros(idstart_loc_t(end)-1,1);
        for i = 1:length(ibins)
            id_targs(idstart_loc_t(i):(idstart_loc_t(i+1)-1)) = idstart_t(ibins(i)):(idstart_t(ibins(i)+1)-1);
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
            ipt_loc = idstart_loc_t(ibinloc):(idstart_loc_t(ibinloc+1)-1);
            id_loc_targ = idstart_t(ibin):(idstart_t(ibin+1)-1);

            
            id_grid_y = j*rbin + (-((npxy-1)/2):((npxy-1)/2))+1+nloc_grid;
            id_grid_loc = (id_grid_x(:).'-1)*(2*nloc_grid+1) + id_grid_y(:);

            sig2targ = pts2grid_t(:,id_loc_targ).' * sig2loc(id_grid_loc(:),:);
                        
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
Asubtract = sparse(iid, jid, vals,size(targ_sort.r(:,:),2),size(src_sort.r(:,:),2));

% toc;
tsub_pre = toc(t1);


end

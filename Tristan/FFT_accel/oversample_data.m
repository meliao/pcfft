function [objout,dataout] = oversample_data(obj,data,nover)
%
%  Oversamples the given geometry given a current discretization
%  and a patchwise oversampling parameter, or a common oversampling
%  order for all patches
%
    npatches = obj.npatches;
    npts = obj.npts;
    if(length(nover)==1)
        noveruse = nover*ones(npatches,1);
    elseif(length(nover)==npatches)
        noveruse = nover;
    else
        fprintf('Incorrect size of npatches, returning same obj\n');
        objout = obj;
        varargout{1} = speye(npts);
    end
    
    ntmp = zeros(npatches,3);
    ntmp(:,1) = obj.norders;
    ntmp(:,2) = noveruse;
    ntmp(:,3) = obj.iptype;
    
    npts_per_patch = zeros(size(noveruse));
    
    npts_per_patch(obj.iptype==1) = ...
       (noveruse(obj.iptype==1)+1).*(noveruse(obj.iptype==1)+2)/2;
   
    npts_per_patch(obj.iptype==11) = ...
       (noveruse(obj.iptype==11)+1).*(noveruse(obj.iptype==11)+1);
   
    npts_per_patch(obj.iptype==12) = ...
       (noveruse(obj.iptype==12)+1).*(noveruse(obj.iptype==12)+1);
   
    npts_per_patch = npts_per_patch(:);
    npts_per_patch = [1; npts_per_patch];
    ixyzso = cumsum(npts_per_patch);
    nptso = ixyzso(end)-1;
    srcover = zeros(12,nptso);
    
    
    [ntmp_uni,~,intmp] = unique(ntmp,'rows');
    nuni = size(ntmp_uni,1);
    
    vmats = cell(nuni,1);
    umats = cell(nuni,1);
    ispl = cell(nuni,1);
    for i=1:nuni
        
        norder = ntmp_uni(i,1);
        nover = ntmp_uni(i,2);
        iptype = ntmp_uni(i,3);
        if(iptype == 1)
            rnodes = koorn.rv_nodes(norder);
            rnodes_over = koorn.rv_nodes(nover);

            vmats{i} = koorn.coefs2vals(norder,rnodes_over);
            umats{i} = koorn.vals2coefs(norder,rnodes);  
            
                
            
        elseif(iptype == 11)
            rnodes = polytens.lege_nodes(norder);
            rnodes_over = polytens.lege_nodes(nover);
            vmats{i} = polytens.lege_coefs2vals(norder,rnodes_over);
            umats{i} = polytens.lege_vals2coefs(norder,rnodes);
        elseif(iptype == 12)
            rnodes = polytens.cheb_nodes(norder);
            rnodes_over = polytens.cheb_nodes(nover);
            vmats{i} = polytens.cheb_coefs2vals(norder,rnodes_over);
            umats{i} = polytens.cheb_vals2coefs(norder,rnodes);
        end
        if(norder>nover)
            [lia,iindb] = ismembertol(rnodes_over',rnodes',1e-7,'ByRows',true);
            ispl{i} = [iindb(lia) find(lia)];
        else
            [lia,iindb] = ismembertol(rnodes',rnodes_over',1e-7,'ByRows',true);
            ispl{i} = [find(lia) iindb(lia)];
        end
    end
       
    dataout = zeros(size(data,1),npatches*(nover+1)*(nover+2)/2);
    
    for i=1:npatches
        istart = ixyzso(i);
        iend = ixyzso(i+1)-1;
        iind = istart:iend;
        srcover(1:9,iind) = obj.srccoefs{i}*vmats{intmp(i)}';
        ru = srcover(4:6,iind);
        rv = srcover(7:9,iind);    
        rtmp = cross(ru,rv);
        vnorm = repmat(vecnorm(rtmp,2),[3,1]);
        srcover(10:12,iind) = rtmp./vnorm;

        istart = obj.ixyzs(i);
        iend = obj.ixyzs(i+1)-1;
        iind2 = istart:iend;
        dataout(:,iind) = data(:,iind2)*umats{intmp(i)}'*vmats{intmp(i)}';

        if(~isempty(ispl{intmp(i)}))
            i1 = ispl{intmp(i)}(:,1);
            i2 = ispl{intmp(i)}(:,2);
            srcover(:,ixyzso(i)+i2-1) = obj.srcvals{i}(:,i1);
        end
    end
    
    objout = surfer(npatches,noveruse,srcover,obj.iptype);
    
    
    
end

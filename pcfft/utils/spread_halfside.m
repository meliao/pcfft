function half_side = spread_halfside(rs, n_nbr, crad)
% hacky way of getting a good radius for the spreading grids
% tries to guarantee that there there is not more than n_nbr*npt
% interactions between points within radius 2*crad*half_side

% 
% % bin points
% T = hypoct(rs,2*n_nbr);
% 
% % for each bin at the finest level
% L = T.nlvl;
% boxes = T.nodes(T.lvp(L)+1:T.lvp(L+1));
% 
% half_sides = zeros(1,length(boxes));
% for i = 1:length(boxes)
%     b = boxes(i);
%     r = 0;
%     for k = 1:size(rs,1)
%         r = r + (rs(k,b.xi) - rs(k,b.xi).').^2;
%     end
%     r = sqrt(r);
% 
%     r = sort(r(:));
% 
%     idr = min(ceil(n_nbr * length(b.xi)), length(r));
% 
%     % get estimated radius as the nth pairwise distance
%     half_side_b = r(idr)/2/crad;
%     half_sides(i) = half_side_b;
% end
% half_side = mean(half_sides);

dim = size(rs,1);
T = hypoct(rs,n_nbr/2^dim);
r = 0;
ileaf = 0;
for i = 1:T.nlvl
    for j = T.lvp(i)+1:T.lvp(i+1)
        if isempty(T.nodes(j).chld)
            r = r + (prod(T.l(:,i))).^(1/dim);
            ileaf = ileaf+1;
        end
    end
end
r = r / ileaf;
r = 2^(1/dim)*2*r;
half_side = r/2/crad;

end
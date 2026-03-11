function half_side = spread_halfside(rs, n_nbr, crad)
% hacky way of getting a good radius for the spreading grids
% tries to guarantee that there there is not more than n_nbr*npt
% interactions between points within radius 2*crad*half_side

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
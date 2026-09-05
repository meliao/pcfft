function half_side = spread_halfside(rs, n_nbr, crad)
% chooses spreading-bin half_side so each point has ~n_nbr points within radius 2*crad*half_side

dim = size(rs,1);
try
    T = hypoct(rs,n_nbr/2^dim);
catch
    error('Error: FLAM not detected. Please install FLAM or manually specify spreading bin halfside (not recommended).')
end
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
r = 2*r/ileaf;
occ_achieved = size(rs,2)/ileaf;
C = gamma(dim/2 + 1)^(1/dim) / (2*sqrt(pi));
half_side = C*r*(n_nbr/occ_achieved)^(1/dim)/(2*crad);

end

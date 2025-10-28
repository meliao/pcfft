function half_side = spread_halfside(rs, n_nbr, crad)
% hacky way of getting a good radius for the spreading grids
% tries to guarantee that there there is not more than n_nbr*npt
% interactions between points within radius 2*crad*half_side


% bin points
T = hypoct(rs,10*n_nbr);
 
% for each bin at the finest level
L = T.nlvl;
boxes = T.nodes(T.lvp(L)+1:T.lvp(L+1));

half_sides = zeros(1,length(boxes));
for i = 1:length(boxes)
    b = boxes(i);
    rx = rs(1,b.xi) - rs(1,b.xi).';
    ry = rs(2,b.xi) - rs(2,b.xi).';
    r = sqrt(rx.^2 + ry.^2);

    r = sort(r(:));
    
    idr = min(ceil(n_nbr * length(b.xi)), length(r));

    % get estimated radius as the nth pairwise distance
    half_side_b = r(idr)/2/crad;
    half_sides(i) = half_side_b;
end
half_side = mean(half_sides);
end
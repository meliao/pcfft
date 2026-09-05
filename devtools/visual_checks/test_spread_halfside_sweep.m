addpath(genpath('../../pcfft'));

nsrc = 5e4;
L = 15;
crad = 2;
n_nbr_list = round(logspace(1, 3, 10));

for dim = [2, 3]
    xs = L*(2*rand(dim,nsrc) - 1);

    tree = createns(xs.');

    actual_nbrs = zeros(size(n_nbr_list));
    half_sides = zeros(size(n_nbr_list));
    for k = 1:numel(n_nbr_list)
        n_nbr = n_nbr_list(k);
        occ = n_nbr/2^dim;
        T = hypoct(xs, occ);
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
        mean_side = r/ileaf;
        half_side = spread_halfside(xs, n_nbr, crad);
        half_sides(k) = half_side;
        radius = 2*crad*half_side;

        idx = randperm(nsrc, min(2000,nsrc));
        counts = rangesearch(tree, xs(:,idx).', radius);
        nc = cellfun(@numel, counts) - 1;
        actual_nbrs(k) = mean(nc);
        fprintf('dim=%d n_nbr=%d occ=%.2f nlvl=%d ileaf=%d nsrc_per_leaf=%.2f mean_leaf_side=%.4f half_side=%.4f radius=%.4f actual=%.2f\n', ...
            dim, n_nbr, occ, T.nlvl, ileaf, nsrc/ileaf, mean_side, half_side, radius, actual_nbrs(k));
    end

    figure(dim); clf
    loglog(n_nbr_list, actual_nbrs, 'o-', 'linewidth', 2)
    hold on
    loglog(n_nbr_list, n_nbr_list, 'k--')
    % loglog(n_nbr_list, actual_nbrs*2^dim, 'o-', 'linewidth', 2)
    hold off
    xlabel('n_{nbr}')
    ylabel('actual mean neighbors')
    title(sprintf('dim = %d', dim))
    legend('measured', 'y = x', 'location', 'northwest')
    grid on
end

function u = apply_eval_over(S,sigma,spread_info, Aquad)
% function u = apply_eval_over(S,sigma,wts,Aspread_s,Aspread_t,Asubtract,kern_hat,isort,isort_t,ngrid, Aquad)

    nover = S.norders(1) + 4;
    [S_over,sigma_over] = oversample_data(S,sigma.',nover);
    wts_over = S_over.wts;

    spread_info.wts = wts_over;
    spread_info.Aquad = sparse(numel(spread_info.isort_t),S_over.npts);

    % u = apply_eval(sigma_over.',wts_over,Aspread_s,Aspread_t,Asubtract,kern_hat,isort,isort_t,ngrid, sparse(numel(isort_t),S_over.npts));
    u = apply_eval(sigma_over.',spread_info);
    u = u + Aquad*sigma;

end
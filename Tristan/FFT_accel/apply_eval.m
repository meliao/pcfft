function u = apply_eval(sigma,spread_info)
% function u = apply_eval(sigma,wts,Aspread_s,Aspread_t,Asubtract,kern_hat,isort,isort_t,ngrid, Aquad)
    sigmasort = sigma(spread_info.isort).*spread_info.wts(spread_info.isort);
    str = spread_info.Aspread_s*sigmasort; str = full(str);

    str_hat = fft2(reshape(str,spread_info.ngrid,spread_info.ngrid));
    u_hat = spread_info.kern_hat.*str_hat;
    ugrid = ifft2(u_hat);
    
    ufast = spread_info.Aspread_t.'*ugrid(:) + spread_info.Asubtract*sigmasort;
    
    uquad = spread_info.Aquad*sigma;
    % ufast = ufast + uquad(isort);

    uquad(spread_info.isort_t) = uquad(spread_info.isort_t) + ufast;
    u = uquad;
    % u = u_hat;
end

function u = apply_mat(sigma,wts,Aspread,Asubtract,kern_hat,isort,ngrid, Aquad)
    sigmasort = sigma(isort).*wts(isort);
    str = Aspread*sigmasort; str = full(str);

    str_hat = fft2(reshape(str,ngrid,ngrid));
    u_hat = kern_hat.*str_hat;
    ugrid = ifft2(u_hat);
    
    ufast = Aspread.'*ugrid(:) + Asubtract*sigmasort;
    
    uquad = Aquad*sigma;
    % ufast = ufast + uquad(isort);

    uquad(isort) = uquad(isort) + ufast;
    u = uquad;
    % u = u_hat;
end


% function u = apply_mat(sigma,Aspread,Asubtract,kern_hat,isort,ngrid)
% tic;
% str = Aspread*sigma(isort); str = full(str);
% tspread = toc
% 
% tic;
% str_hat = fft2(reshape(str,ngrid,ngrid));
% u_hat = kern_hat.*str_hat;
% ugrid = ifft2(u_hat);
% tconv = toc
% 
% tic;
% u = Aspread.'*ugrid(:);
% tinterp = toc
% 
% tic;
% u = u + Asubtract*sigma(isort);
% tadd_sub = toc
% end

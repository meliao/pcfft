% % N = 1000; % number of sources
% % M = 3000; % number of targets
% % src_info = struct();
% % targ_info = struct();
% % src_info.r = rand(3, N);
% % src_info.kappa = rand(1, N);
% % targ_info.r = rand(3, M);
% % targ_info.n = rand(3, M);
% % 
% % a = kern(src_info,targ_info);
% % size(a)
% % 
% % function evals = kern(src_info, targ_info)
% %     % evaluate the gradient of the electrostatic kernel between N source points src_info.r and
% %     % M target points targ_info.r
% %     %
% %     % Output shape is (3M, N)
% % 
% %     % Shape (M, N)
% %     rx = src_info.r(1, :) - targ_info.r(1, :).';
% %     ry = src_info.r(2, :) - targ_info.r(2, :).';
% %     rz = src_info.r(3, :) - targ_info.r(3, :).';
% % 
% %     dist = sqrt(rx.^2 + ry.^2 + rz.^2);
% %     gradx = rx ./ dist.^2; gradx = reshape(gradx,1,[]);
% %     grady = ry ./ dist.^2; grady = reshape(grady,1,[]);
% %     gradz = rz ./ dist.^2; gradz = reshape(gradz,1,[]);
% % 
% %     % Shape (3M, N)
% %     evals = reshape([gradx;grady;gradz], 3*size(targ_info.r(:,:),2), size(src_info.r(:,:),2));
% %     evals = evals / (4*pi);
% % end
% 
% 
% 
% N = 1000; % number of sources
% M = 3000; % number of targets
% src_info = struct();
% src_info.r = rand(2, N);
% src_info.n = randn(2, N); % normal vector at each source point
% targ_info = struct();
% targ_info.r = rand(2, M);
% 
% 
% tol = 1e-6;
% [grid_info, proxy_info] = get_grid(@kern, src_info, targ_info, tol);
% 
% [A_spread_src, srt_info_src] = get_spread(@kern, @kern_s, src_info, ...
%                                           grid_info, proxy_info, {'r','n'});
% [A_spread_targ, srt_info_targ] = get_spread(@kern, [], targ_info, ...
%                                              grid_info, proxy_info);
% 
% function k_evals = kern(src_pts, target_pts)
%    % src_pts is a struct where src_pts.r with shape (2, M)
%    % target_pts is a struct where target_pts.r has shape (2, N)
%    % Computes log{|| src - target||}
%    % Output shape is (N, M)
% 
%    % Shape (N, M)
%    rx = src_pts.r(1, :) - target_pts.r(1, :).';
%    ry = src_pts.r(2, :) - target_pts.r(2, :).';
% 
%    dist = sqrt(rx.^2 + ry.^2);
% 
%    k_evals = log(dist);
% end
% 
% function k_evals = kern_s(src_pts, target_pts)
%    % src_pts is a struct where
%    %         - src_pts.r with shape (2, M)
%    %         - src_pts.n with shape (2, M)
%    % target_pts is a struct where  target_pts.r has shape (2, N)
%    % Computes \partial_{n(src)}log{|| src - target||}
%    % Output shape is (N, M)
% 
%    % Shape (N, M)
%    rx = src_pts.r(1, :) - target_pts.r(1, :).';
%    ry = src_pts.r(2, :) - target_pts.r(2, :).';
% 
%    dist = sqrt(rx.^2 + ry.^2);
% 
%    k_evals = src_pts.n(1,:).*rx + src_pts.n(2,:).*ry;
%    k_evals = k_evals ./ dist.^2;
% end

% Arsc = rscMat(A_addsub_c);
b = randn(size(A_addsub_c,1),1);
tic;c = A_addsub_c*b; toc;


% tic;d = Arsc*b; toc;

tic; sys_app(b);toc;
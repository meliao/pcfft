function my_scatter(x,varargin)
% if isempty(varargin)
%     scatter3(x(1,:), x(2,:), x(3,:))
% else
%     scatter3(x(1,:), x(2,:), x(3,:),varargin{:})
% end

if isempty(varargin)
    scatter(x(1,:), x(2,:))
else
    scatter(x(1,:), x(2,:),varargin{:})
end

end

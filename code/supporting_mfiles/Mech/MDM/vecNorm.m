function v  = vecNorm(u,dim)
% v  = vecNorm(u,dim)
% v = u/sqrt(ui^2);

if nargin < 2, dim = 1; end

if size(u,1) == 1, dim = 1; end
if size(u,2) == 2, dim = 2; end

v = bsxfun(@rdivide, u, sqrt(sum(u.^2,dim)));

end








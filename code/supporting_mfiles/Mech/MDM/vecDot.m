function w = vecDot(u,v,dim)
% w = vecDot(u,v,dim)
% w = u.v

if nargin < 3, dim = 2; end
w = sum(u.*v,dim);


end
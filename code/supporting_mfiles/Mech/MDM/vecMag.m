function mag = vecMag(u,dim)
% mag  = vecMag(u,dim)
% mag = sqrt(ui^2)

if nargin < 2, dim = 1; end

if size(u,1) == 1, dim = 1; end
if size(u,2) == 2, dim = 2; end

mag = sqrt(sum(u.^2,dim));

end
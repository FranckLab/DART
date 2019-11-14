function cosTheta = vecCos(u,v,dim)
% cosTheta  = vecNorm(u,v,dim)
% cosTheta = u.v/(|u||v|)
if nargin < 3, dim = 1; end

if size(u,1) == 1 || size(v,1) == 1, dim = 2; end
if size(u,2) == 2 || size(v,2) == 2, dim = 1; end



cosTheta = sum(bsxfun(@times,u,v),dim);

cosTheta = cosTheta./(  sqrt(sum(u.^2,dim)).*sqrt(sum(v.^2,dim))  ); 

end

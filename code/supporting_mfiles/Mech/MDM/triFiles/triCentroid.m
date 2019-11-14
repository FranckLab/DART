function [CA, CV, dCA, dCV] = triCentroid(f,v,w)
% v = vertex coordinates
% f = triangulation
% w = weights for each vertex or for each triangle (if no input then
%     calculate area centroid)

v0 = v(f(:,1),:);
v1 = v(f(:,2),:);
v2 = v(f(:,3),:);
v3 = cross(v1 - v0,v2 - v0);  % normal vector (not normalized)

dA = 0.5*sqrt(sum(v3.^2,2));
A = sum(dA);


dV = 1/6*abs(sum(v3.*(bsxfun(@minus, v0,  mean(v0))),2)); % volume for each pyramid
V = sum(dV);
 
if nargin < 3 % calculate weights based on area
    
    dCA = bsxfun(@times, (v0+v1+v2)/3,dA);
    CA = sum(dCA)/A;
    dCV = bsxfun(@times, (v0 + v1 + v2 + (v3 + v0))/4, dV);
    CV = sum(dCV)/V;
    
elseif size(w,1) == size(v,1) % weight given for each vertex
    
    w0 = w(f(:,1)); % get weight for each vertex
    w1 = w(f(:,2));
    w2 = w(f(:,3));
   
    W = mean(w0 + w1 + w2)/3; % mean weight of triangles
%     W = sum(w0 + w1 + w2)/3;
    dCA = (repmat(w0,[1,3]).*v0+repmat(w1,[1,3]).*v1+repmat(w2,[1,3]).*v2)/3;
    dCA = bsxfun(@times, dCA, dA);
    CA = sum(dCA)/(W*A);
% C = sum(dC)/W;
end

end
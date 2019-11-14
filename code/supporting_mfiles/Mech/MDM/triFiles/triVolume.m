
function [V, dV, C] = triVolume(f,v)
% [V, dV] = triVolume(f,v)
% f = triangulation
% v = vertex coordinates

v = bsxfun(@minus, v, mean(v));
v0 = v(f(:,1),:);
v1 = v(f(:,2),:);
v2 = v(f(:,3),:);

v3 = cross(v1 - v0,v2 - v0);  % normal vector (not normalized)
dV = 1/6*abs(sum(v3.*v0,2)); % volume for each pyramid
V = sum(dV);  % total volume

% C = sum((v0 + v1 + v2 + v3)/4)/V;
end
function [A, dA] = triArea(f,v)
% v = vertex coordinates
% f = triangulation


v0 = v(f(:,1),:);
v1 = v(f(:,2),:);
v2 = v(f(:,3),:);
t1 = v1 - v0;
t2 = v2 - v0;

% calculation triangulation area
dA = 0.5*sqrt(sum(cross(t1,t2).^2,2));
A = sum(dA);
end
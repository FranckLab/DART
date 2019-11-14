function [V,A] = stlVolume(f,v)
% Given a surface triangulation, compute the volume enclosed using
% divergence theorem.
% Assumption:Triangle nodes are ordered correctly, i.e.,computed normal is outwards
% Input: p: (3xnPoints), t: (3xnTriangles)
% Output: total volume enclosed, and total area of surface  
% Author: K. Suresh; suresh@engr.wisc.edu

% Compute the vectors d13 and d12
d13= [(v(1,t(2,:))-p(1,f(3,:))); (p(2,t(2,:))-p(2,t(3,:)));  (p(3,t(2,:))-p(3,t(3,:)))];
d12= [(p(1,t(1,:))-p(1,t(2,:))); (p(2,t(1,:))-p(2,t(2,:))); (p(3,t(1,:))-p(3,t(2,:)))];
cr = cross(d13,d12,1);%cross-product (vectorized)
dA = 0.5*sqrt(cr(1,:).^2+cr(2,:).^2+cr(3,:).^2);% Area of each triangle
A = sum(dA);
crNorm = sqrt(cr(1,:).^2+cr(2,:).^2+cr(3,:).^2);
zMean = (p(3,t(1,:))+p(3,t(2,:))+p(3,t(3,:)))/3;
nz = -cr(3,:)./crNorm;% z component of normal for each triangle
dV = dA.*zMean.*nz; % contribution of each triangle
V = sum(dV);%divergence theorem
%%

% 'z' coordinate of triangle centers. 
FaceCentroidZ = ( V(F(:, 1), 3) + V(F(:, 2), 3) + V(F(:, 3), 3) ) /3; 
% Face normal vectors, with length equal to triangle area. 
FNdA = cross( (V(F(:, 2), :) - V(F(:, 1), :)), ... 
(V(F(:, 3), :) - V(F(:, 2), :)) , 2 ) / 2; 
% Volume from divergence theorem (using vector field along z). 
V = FaceCentroidZ' * FNdA(:, 3);
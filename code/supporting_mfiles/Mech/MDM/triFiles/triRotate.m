function v = triRotate(v,R)
% v = triRotate(v,R);
% R can be a [3x3] rotation matrix or axis-angle representation 
% [ux, uy, uz, angle (radians)]


if numel(R) == 4
    R = vrrotvec2mat(R);
end

C = mean(v);
v = bsxfun(@minus, v, C);
    
 for i = 1:size(v,1)
    v(i,:) = R*v(i,:)';
 end

v = bsxfun(@plus, v, C);

end
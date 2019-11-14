function v = triAffine(v,A)
% v = triScale(v,A);
% A is 3x3 affine transform

C = mean(v);
v = bsxfun(@minus, v, C);

for i = 1:size(v,1)
    v(i,:) = v(i,:)*A;
end

v = bsxfun(@plus, v, C);

end
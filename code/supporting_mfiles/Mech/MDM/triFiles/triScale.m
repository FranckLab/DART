function v = triScale(v,A)
% v = triScale(v,A);
% A can be a [3x3] scaling matrix
% [Ax  0   0 ]
% [0   Ay  0 ]
% [0    0  Az];
% or vector, A = [Ax, Ay, Az]
% or scalar, A = Ax = Ay = Az;

if numel(A) == 3, A = diag(A); end
if numel(A) == 1, A = diag([A A A]); end

if ~all(diag(A) == [1 1 1]')
    
    C = mean(v);
    v = bsxfun(@minus, v, C);
    
    for i = 1:size(v,1)
        v(i,:) = v(i,:)*A;
    end
    
    v = bsxfun(@plus, v, C);
end

end



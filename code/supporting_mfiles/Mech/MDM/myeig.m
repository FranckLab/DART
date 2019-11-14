function [V D]=myeig(X,varargin)
% Find eigen vector/value and enforce a sign convention by
% making the largest eigen vector element non-negative.
[V D]=eig(X,varargin{:});
[~,I] = max(abs(V),[],1);
I = sub2ind(size(V),I,1:size(V,2));
V = bsxfun(@times,V,sign(V(I)));
end % myeig
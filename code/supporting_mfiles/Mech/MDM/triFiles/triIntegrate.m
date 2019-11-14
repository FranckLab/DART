function [int, A, dA] = triIntegrate(f,v,w)
% [int, A, dA] = triIntegrate(f,v,w)
% v = vertex coordinates
% f = triangulation
% w = weights for each vertex or for each triangle (if no input then
%     calculate area centroid)



v0 = v(f(:,1),:);
v1 = v(f(:,2),:);
v2 = v(f(:,3),:);

dA = 0.5*sqrt(sum(cross(v1 - v0,v2 - v0).^2,2));
A = sum(dA);
int = A;

if nargin > 2,
    
    if size(w,1) == size(v,1) % weight given for each vertex
        
        for i = 1:size(w,2)
            w0 = w(f(:,1),i); % get weight for each vertex
            w1 = w(f(:,2),i);
            w2 = w(f(:,3),i);
            int(i) = sum((w0 + w1 + w2)/3.*dA);  % int(w,dA)
        end
        
    elseif size(w,1) == size(f,1)
        
        for i = 1:size(w,2)
            int(i) = sum(w(:,i).*dA);  % int(w,dA)
        end
        
    elseif size(w,3) == size(f,1)
        
        for i = 1:size(w,1)
            for j = 1:size(w,2)
                int(i,j) = sum(reshape(w(i,j,:),[],1).*dA);
            end
        end
        
    end
    
end




end
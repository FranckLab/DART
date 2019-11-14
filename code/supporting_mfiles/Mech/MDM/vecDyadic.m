function W = vecDyadic(u,v)

u = u';

W = zeros(size(u,1),size(u,1),size(u,2));

for i = 1:size(W,3)
   W(:,:,i) = u(:,i)*v(i,:); 
end

end
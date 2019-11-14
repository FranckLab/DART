function [dA, cX, cY, cZ] = dASurface(sX,sY,sZ)

kernalSize = 2;
xWindows = im2col(sX,kernalSize*[1 1],'sliding');
yWindows = im2col(sY,kernalSize*[1 1],'sliding');
zWindows = im2col(sZ,kernalSize*[1 1],'sliding');

dA = zeros(1,size(xWindows,2));
cX = zeros(1,size(xWindows,2));
cY = zeros(1,size(xWindows,2));
cZ = zeros(1,size(xWindows,2));

for i = 1:size(xWindows,2)
    x_ = xWindows(:,i);
    y_ = yWindows(:,i);
    z_ = zWindows(:,i);
    
    [dA(i), cX(i), cY(i), cZ(i)] = trapArea(x_,y_,z_);
end

dA = reshape(dA,(size(sX)-1));
cX = reshape(cX,(size(sX)-1));
cY = reshape(cY,(size(sX)-1));
cZ = reshape(cZ,(size(sX)-1));
end

function [dA, Cx, Cy, Cz] = trapArea(x_,y_,z_)

dt = DelaunayTri(x_,y_);
xyz = [x_, y_, z_];
triRepTrap = TriRep(dt.Triangulation,xyz);
tri = dt.Triangulation;
C = incenters(triRepTrap);

dA1 = triArea(tri(1,:),xyz);
dA2 = triArea(tri(2,:),xyz);

dA = dA1 + dA2;

Cx = (C(1,1)*dA1 + C(2,1)*dA2)/dA;
Cy = (C(1,2)*dA1 + C(2,2)*dA2)/dA;
Cz = (C(1,3)*dA1 + C(2,3)*dA2)/dA;


trimesh(triRepTrap)
hold on;
plot3(C(1,1),C(1,2),C(1,3),'ro');
plot3(C(2,1),C(2,2),C(2,3),'*');
plot3(Cx,Cy,Cz,'o')
box on;
hold off


end

function dA = triArea(tri,xyz)
xyz0 = xyz(tri(1),:);
xyz1 = xyz(tri(2),:);
xyz2 = xyz(tri(3),:);

t1 = xyz1 - xyz0;
t2 = xyz2 - xyz0;

dA = 0.5*norm(cross(t1,t2));
end
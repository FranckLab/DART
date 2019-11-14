function [dA, dtTri, centers] = dVTriangulation(s)

sX_ = s{1}(:); sY_ = s{2}(:); sZ_ = s{3}(:);
dt = DelaunayTri([sX_,sY_,sZ_]);

dtTri = dt.Triangulation;
s = [sX_, sY_, sZ_];
triS = TriRep(dtTri,s);
centers = incenters(triS);


xyz0 = s(dtTri(:,1),:);
xyz1 = s(dtTri(:,2),:);
xyz2 = s(dtTri(:,3),:);
t1 = xyz1 - xyz0;
t2 = xyz2 - xyz0;
dA = 0.5*sqrt(sum(cross(t1,t2).^2,2));


end
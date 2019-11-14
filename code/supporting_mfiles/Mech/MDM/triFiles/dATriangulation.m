function tri = dATriangulation(s)

sX_ = s{1}(:); sY_ = s{2}(:); sZ_ = s{3}(:);
dt = delaunayTriangulation(sX_,sY_);

dtTri = dt.ConnectivityList;
s_ = [sX_, sY_, sZ_];
triS = triangulation(dtTri,s_);
centers = incenter(triS);


xyz0 = s_(dtTri(:,1),:);
xyz1 = s_(dtTri(:,2),:);
xyz2 = s_(dtTri(:,3),:);
t1 = xyz1 - xyz0;
t2 = xyz2 - xyz0;
dA = 0.5*sqrt(sum(cross(t1,t2).^2,2));

tri.dA = dA;
tri.triangulation = dtTri;
tri.centers = centers;
tri.vertices = dt.Points;
end
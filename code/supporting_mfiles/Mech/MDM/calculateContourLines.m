function [polys, planesOrigin] = calculateContourLines(f,v,m, nPlanes)
%
% m = central axis vector

C = triCentroid(f,v);
m = m/norm(m);

line = [C m];

[intP, intP(:,4)] = intersectLineMesh3d(line, v, f);
intP = [intP(1,:); intP(end,:)];

planesOrigin = zeros(nPlanes,3);
for i = 1:4
    planesOrigin(:,i) = linspace(intP(1,i),intP(2,i), nPlanes);
end

planes = zeros(nPlanes,9);
for i = 1:nPlanes
    planes(i,:) = createPlane(planesOrigin(i,1:3),m);
    polys(i,:) = intersectPlaneMesh(planes(i,:), v, f);
end

% 
% 
% figure;
% hold on;
% p = patch('Faces',f,'Vertices',v,'FaceColor','c');
% set(p, 'EdgeAlpha', 0.10,'facealpha',0.50);
% % trimesh(f,v(:,1),v(:,2),v(:,3),'facealpha',0.50,'edgealpha',0.10, 'FaceColor', 'c');
% plot3(C(1),C(2),C(3),'*','linewidth',3);
% plot3(intP(:,1),intP(:,2),intP(:,3),'ro');
% drawLine3d(line,'color','r');
% % for i = 1:nPlanes, drawPlane3d(planes(i,:)); end
% for i = 1:nPlanes, drawPolygon3d(polys(i,:)); end
% hold off;
% axis image; 
% view(3)
% light
% set(gca,'color','none','xtick',[],'ytick',[],'ztick',[])
%     box on;
%

end
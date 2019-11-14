function [fE, vE, ellipseInfo] = triMinEllipse(v, edgeLength, plot01,varargin)

if nargin < 3, plot01 = false; end
if nargin < 2, edgeLength = 0.25; end

fCH = convhull(v);

% remove points within CH
v0CH = v(unique(fCH(:)),:);
[A, C] = MinVolEllipse(v0CH', .01);
[~, R, rotation] = svd(A);
R = 1./sqrt(diag(R))';

ellipseFun =@(p) (p(:,1).^2/R(1).^2+p(:,2).^2/R(2).^2+p(:,3).^2/R(3)^2-1);
[vE, fE]=distmeshsurface(ellipseFun,@huniform,edgeLength,1.1*[-R;R],10);

vE = bsxfun(@plus, vE, C');
vE = triRotate(vE,rotation);

ellipseInfo.centerForm = A;
ellipseInfo.radius = R;
ellipseInfo.rotation = rotation;
ellipseInfo.centroid = C;
ellipseInfo.function = ellipseFun;


if plot01
    figure;
    hold on;
%     trisurf(fCH,v(:,1),v(:,2),v(:,3) ,  'FaceColor', 'b');
    trisurf(fE,vE(:,1),vE(:,2),vE(:,3),   'FaceColor', 'c', 'faceAlpha', 0.2,'edgealpha',0.05);
    light
%         plot3(v0CH(:,1),v0CH(:,2),v0CH(:,3),'.k')
    trisurf(f,v(:,1),v(:,2),v(:,3),   'FaceColor', 'g', 'edgecolor','w', 'faceAlpha', 1,'edgealpha',0.2);
    lighting phong
    view(3);
    set(gca,'color','none','xtick',[],'ytick',[],'ztick',[])
    box on;
    axis image
    hold off;
    
%     export_fig(gcf,'I.png','-png','-nocrop','-opengl','-r200');
end

end


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

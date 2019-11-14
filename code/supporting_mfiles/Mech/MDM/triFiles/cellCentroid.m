function [C, CEp] = cellCentroid(f,v,Ep,dm,plot01)
% f = faces;
% v = vertices
% Ep = principal strains defined for each face
% dm = meshspacing
% plot01 = plot boolean

%% LOAD DATA
% clear; clc;
% % load('Celldata_dt02.mat','FV','vertices');
% load('psdt2.mat');
% Ep{1} = psdt2E1;
% Ep{2} = psdt2E2;
% Ep{3} = psdt2E3;
% f = FV.faces;
% v = vertices;
% clearvars -except Ep f v


%% INTERPOLATE PRINCIPAL STRAINS ON SURFACE
% create meshgrid
mSize = size(Ep{1});
mIdx = cell(1,3);
for i = 1:3, mIdx{i} = 1:dm:dm*mSize(i); end
[m{1}, m{2}, m{3}] = ndgrid(mIdx{:});
EpV = zeros(size(v));

for i = 1:3
    F = griddedInterpolant(m{1},m{2},m{3},Ep{i},'linear');
    EpV(:,i) = F(v(:,2),v(:,1),v(:,3));
end
EpV(:,4) = sqrt(sum(EpV.^2,2));

[CEp, ~] = triCentroid(v,f,EpV(:,4));


%% GET CENTROIDS
[C, ~] = triCentroid(v,f); % get volume centroid

Cvec = CEp - C;
% velocity = [10 0 0];

if plot01
    figure;
    hold on;
    % h = drawMesh(v,f);
    h = trisurf(f,v(:,1),v(:,2),v(:,3),EpV(:,4));
    set(h,'linestyle','none')
    alpha(0.25)
    view(3); axis('vis3d');
    title('Icosahedron');
    plot3(C(1),C(2),C(3),'o','MarkerSize',10,'MarkerEdgeColor','g','MarkerFaceColor','g')
    plot3(CEp(1),CEp(2),CEp(3),'o','MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor','r')
    quiver3(C(1),C(2),C(3),Cvec(1),Cvec(2),Cvec(3),1,'linewidth',2)
    % quiver3(C(1),C(2),C(3),velocity(1),velocity(2),velocity(3),1,'linewidth',2)
    lighting gouraud
    hold off;
end

end


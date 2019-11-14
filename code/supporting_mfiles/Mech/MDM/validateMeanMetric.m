%% runValidated Neutrophil 3D
close all; clear

sphR = 32;
mSize = [256,256,256];

[f, v] = sphere_tri('ico',4,sphR);
triC = [0 0 0];

v = bsxfun(@plus, v, triC + (mSize-1)/2);
[n,~,p] = triNormals(f,v);
n = -1*n;

triV = triVolume(f, v);

% for i = 1:3, mIdx{i} = (1:mSize(i)) - (mSize(i)+1)/2; end
for i = 1:3, mIdx{i} = (1:mSize(i)); end
[m{1}, m{2}, m{3}] = meshgrid(mIdx{:});
%%
% uniaxial stretch
% L(3) = 1.1;
% L(1) = 1/sqrt(L(3));
% L(2) = 1/sqrt(L(3));
% 
% u{1} = (L(1) - 1)*m{1};
% u{2} = (L(2) - 1)*m{2};
% u{3} = (L(3) - 1)*m{3};

% shear
% k = 0.1*tand(45);
%
% u{1} = k*m{2};
% u{2} = zeros(mSize);
% u{3} = zeros(mSize);

% biaxial stretch
% L(1) = 1.1;
% L(2) = 0.8;
% L(3) = 1/(L(1)*L(2));
% 
% u{1} = (L(1) - 1)*m{1};
% u{2} = (L(2) - 1)*m{2};
% u{3} = (L(3) - 1)*m{3};

% translation with strech
% T = [1 3 10];
% L(3) = 1.1;
% L(1) = 1/sqrt(L(3));
% L(2) = 1/sqrt(L(3));
% 
% u{1} = (L(1) - 1)*(m{1} - (mSize(1)-1)/2 - triC(1));
% u{2} = (L(2) - 1)*(m{2} - (mSize(1)-1)/2 - triC(2));
% u{3} = (L(3) - 1)*(m{3} - (mSize(1)-1)/2 - triC(3));

%Pure rotation around an axis
RAxis = [1 0 0];
RAngle = 45*pi/180;
R =  vrrotvec2mat([RAxis RAngle]);

I = eye(3,3);
%u = (Rij-dij) xj
for i=3:-1:1
   u{i} = (R(i,1)-I(i,1))*(m{1} - (mSize(1)-1)/2 - triC(1)) + (R(i,2)-I(i,2))*(m{2} - (mSize(1)-1)/2 - triC(2)) +(R(i,3)-I(i,3))*(m{3} - (mSize(1)-1)/2 - triC(3));
end

uTri = triInterp(p,u,1);
%uTri = bsxfun(@plus, uTri, T);
vDef = triInterp(v,u,1) + v;
%vDef = bsxfun(@plus, vDef, T);


[~, ~, cDef] = triNormals(f,vDef);
triVDef = triVolume(f,vDef);

triV = triVolume(f,v);

ui_nj = vecDyadic(uTri,n);
F = I + (1/triV)*triIntegrate(f,v, ui_nj)
J = det(F);
C = F'*F;
B = F*F';
U = sqrtm(C);
V = sqrtm(B);

E = 1/2*(C - I);
e = triIntegrate(f,v, ui_nj);
e = 1/(2*triV)*(e + e');
RMean = V\F;
axisAngleR = vrrotmat2vec(RMean);
[eigVecE, eigValE] = eig(E);
eigValE = diag(eigValE);
% 
% clc
% x = p;
% xi_vj = vecDyadic(x,uTri);
%
% for j = 1:size(xi_vj,3)
%     xi_vj_nj(j,:) = xi_vj(:,:,j)*n(j,:)';
% end
%
% vCM = 1/triV*triIntegrate(f,v, xi_vj_nj)



% x = bsxfun(@minus, p, triCentroid(f,v));
% xi_vj = vecDyadic(x,uTri);
% 
% for j = 1:size(xi_vj,3)
%     xi_vj_nj(j,:) = xi_vj(:,:,j)*n(j,:)';
% end
% vCM = 1/triV*triIntegrate(f,v, xi_vj_nj)

%% u in the volume
figure;

view_ = {2, [90, 0], 3};

%for i = 1:3
%    subplot(1,3,i)
hold on;
h{1} = trisurf(f,v(:,1),v(:,2),v(:,3));
%h{2} = trisurf(f,vDef(:,1),vDef(:,2),vDef(:,3));
% plot3(c(:,1),c(:,2),c(:,3),'.');
unorm = sqrt(u{1}.^2+u{2}.^2+u{3}.^2);
maxunorm = max(unorm(:));

skip = 16;
h{2} = quiver3(m{1}(1:skip:end,1:skip:end,1:skip:end),m{2}(1:skip:end,1:skip:end,1:skip:end),m{3}(1:skip:end,1:skip:end,1:skip:end),   u{1}(1:skip:end,1:skip:end,1:skip:end), u{2}(1:skip:end,1:skip:end,1:skip:end), u{3}(1:skip:end,1:skip:end,1:skip:end), 1);
hold off

axis equal
%axis(1.25*[-sphereRadius sphereRadius -sphereRadius sphereRadius -sphereRadius sphereRadius]);
xlabel('x');
ylabel('y');
zlabel('z');
% colorbar
%view(view_{i})
% set(h{1},'linestyle','none');
%set(h{2},'linestyle','none','edgecolor', 0.90*[1 1 1]);
%set(h{1}, 'FaceColor', 'w')
%set(h{2}, 'FaceColor', 'w')
%alpha(h{2},0.25)
set(gca,'Color',[0 0 0]);
% box(gca)

%%
figure;

tit{1} = ['Applied   \lambda: [',num2str(diag(F)','%1.3f '),']'];
tit{2} = ['Recovered \lambda: [',num2str(L,'%1.3f '),']' ];
tit{3} = '------------------------------------------------------';
tit{4} = ['Applied   <v_c_m> = [', num2str(T,'%1.1f '),']'];
tit{5} = ['Recovered <v_c_m> = [', num2str(vCM,'%1.1f '),']'];

view_ = {2, [90, 0], 3};

for i = 1:3
    subplot(1,3,i)
    hold on;
%     h{1} = trisurf(f,v(:,1),v(:,2),v(:,3));
    h{2} = trisurf(f,vDef(:,1),vDef(:,2),vDef(:,3));
    % plot3(c(:,1),c(:,2),c(:,3),'.');
%     quiver3(p(:,1),p(:,2),p(:,3),   uTri(:,1),uTri(:,2),uTri(:,3), 0);
    hold off
    
    
    axis equal
    axis(reshape([mean(v) - 1.5*sphR ;mean(v) + 1.5*sphR],1,[]));
    xlabel('x');
    ylabel('y');
    zlabel('z');
    % colorbar
    view(view_{i})
    % set(h{1},'linestyle','none');
%     set(h{1}, 'linestyle', 'none')
    set(h{2}, 'linestyle','none');
    
%     set(h{1}, 'FaceColor', 'w')
%     set(h{2}, 'FaceColor', 'w')
%     alpha(h{2},0.25)
    set(gca,'Color',[0 0 0]);
    
    light
    % box(gca)
    
    if i == 2,
        title(tit,'fontsize',16)
    end
    
    
end




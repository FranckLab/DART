%% runValidated Neutrophil 3D
close all; clear

sphereRadius = 32;
mSize = 64*[1 1 1];

[f, v] = sphere_tri('ico',4,sphereRadius);

sc=1;
v = bsxfun(@times, v, [1,sc,sc]);
[fn,~,fp] = triNormals(f,v);

triV = triVolume(f, v);

for i = 1:3, mIdx{i} = (1:mSize(i)) - (mSize(i)+1)/2; end
[m{1}, m{2}, m{3}] = meshgrid(mIdx{:});


% uniaxial stretch
% L{3} = [linspace(1.1,1.0,mSize(3)/2), linspace(1.0,1.1,mSize(3)/2)];
%%
% % uniaxial stretch
L(3) = 1.2;
L(2) = 1/sqrt(L(3));
L(1) = 1/sqrt(L(3));

u{1} = (L(1) - 1)*m{1};
u{2} = (L(2) - 1)*m{2};
u{3} = (L(3) - 1)*m{3};

% L(3) = 1.1;
% L(2) = 1;
% L(1) = 1;
% 
% u{1} = (L(1) - 1)*m{1};
% u{2} = (L(2) - 1)*m{2};
% u{3} = (L(3) - 1)*m{3};

% % shear
% k = 0.1*tand(45);
% 
% u{1} = k*m{2};
% u{2} = zeros(mSize);
% u{3} = zeros(mSize);
%
% % biaxial stretch
% L(1) = 1.1;
% L(2) = 0.8;
% L(3) = 1/(L(1)*L(2));
% %
% u{1} = (L(1) - 1)*m{1};
% u{2} = (L(2) - 1)*m{2};
% u{3} = (L(3) - 1)*m{3};

dm = 1;
Fij = calculateFij({u},dm,'prewitt');
FTri = triInterp(fp,Fij{1},mIdx);
ETri = calculateEijTri(FTri);

for i = 1:size(ETri,3)
    [eigVec_, eigVal_] = eig(ETri(:,:,i));
    eigVal_ = diag(eigVal_)';
    minEigVal(i,:) = [eigVec_(:,1)', eigVal_(1)];
    maxEigVal(i,:) = [eigVec_(:,3)', eigVal_(3)];
    
    eigValNorm(i,:) = eigVal_/norm(eigVal_);
    
    
    eigVal(i,:) = [eigVec_(:,1)', eigVal_(1), eigVec_(:,3)', eigVal_(3)];
    %             eigValW(ii,:) = [eigVec_(:,1)'*eigVal_(1,1), eigVal_(1,1), eigVec_(:,3)'*eigVal_(3,3), eigVal_(3,3)];
    
end

% threshold min
min01 = eigValNorm(:,1) > prctile(eigValNorm(:,1),25) & eigValNorm(:,1) < prctile(eigValNorm(:,1),75);
max01 = eigValNorm(:,3) > prctile(eigValNorm(:,3),25) & eigValNorm(:,3) < prctile(eigValNorm(:,3),75);
%%
dataInfo.varnames = {'x1','x2','x3','vmin1','vmin2','vmin3','Emin', 'vmax1','vmax2','vmax3','Emax'};
dataInfo.zonename = 'Cell_';
TPData = convert2Tecplot('FESurfaces', [{f},{v}],  {eigVal}, dataInfo);
mat2tecplot(TPData,'C:\Users\Eyal\Dropbox\Neutrophil3D\mFiles\E.plt');
%% 


stream3(x,y,z,u,v,w,sx,sy,sz)




%%


% eigValMag = vecMag(eigValETri,2);

figure;
view_ = {2, [90, 0], 3};
for i = 1:3
    
    subplot(1,3,i)
%     h{1} = trisurf(f,v(:,1),v(:,2),v(:,3),maxEigVal(:,4));
    h{1} = trisurf(f,v(:,1),v(:,2),v(:,3),reshape(ETri(3,3,:),1,[]));
    
    
    axis(1.25*[-sphereRadius sphereRadius -sphereRadius sphereRadius -sphereRadius sphereRadius]);
    xlabel('x');
    
    ylabel('y');
    zlabel('z');
    
    view(view_{i})
    set(h{1},'linestyle','none');
    set(gca,'Color',[0 0 0]);
    % box(gca)
    
    light
    axis equal
    
    if i == 2, colorbar('location','south'); end
    
end


%%



plotRange = 1.1*[min(minEig) max(maxEig)];


figure
hold on;
plot(bins, maxEig, 'r-','linewidth', 3)
plot(bins, minEig, 'b-','linewidth', 3)
plot([0 0], [-1 1], 'k--', 'linewidth', 2);
hold off

axis([-1,1,plotRange(1) plotRange(2)])

fontSize = 10;
xlabel('Normalized Distance along velocity axis','FontSize',fontSize);
ylabel('Principal Strain','FontSize',fontSize);
set(gca,'fontsize',fontSize)
hl = legend('max(Ei)','max(Ei)');
set(hl,'FontSize',fontSize);
box on;
drawnow









%%
figure;

view_ = {2, [90, 0], 3};

for i = 1:3
    subplot(1,3,i)
    hold on;
    h{1} = trisurf(f,v(:,1),v(:,2),v(:,3));
    h{2} = trisurf(f,vDef(:,1),vDef(:,2),vDef(:,3));
    % plot3(c(:,1),c(:,2),c(:,3),'.');
    quiver3(fp(:,1),fp(:,2),fp(:,3),   uTri(:,1),uTri(:,2),uTri(:,3), 0);
    hold off
    
    
    axis equal
    axis(1.25*[-sphereRadius sphereRadius -sphereRadius sphereRadius -sphereRadius sphereRadius]);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    % colorbar
    view(view_{i})
    % set(h{1},'linestyle','none');
    set(h{2},'linestyle','none','edgecolor', 0.90*[1 1 1]);
    set(h{1}, 'FaceColor', 'w')
    set(h{2}, 'FaceColor', 'w')
    alpha(h{2},0.25)
    set(gca,'Color',[0 0 0]);
    % box(gca)
    
    
end
%%





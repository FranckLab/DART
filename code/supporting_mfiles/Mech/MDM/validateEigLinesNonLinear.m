%% runValidated Neutrophil 3D
close all; clear

sphereRadius = 16;
mSize = [128,128,128];

[f, v] = sphere_tri('ico',4,sphereRadius);
[fn,~,fp] = triNormals(f,v);

triV = triVolume(f, v);

for i = 1:3, mIdx{i} = (1:mSize(i)) - (mSize(i)+1)/2 ; end
[m{1}, m{2}, m{3}] = meshgrid(mIdx{:});


%%
% L{3} = [linspace(0.8,1.2,mSize(3))];
% L{1} = 1;
% L{2} = 1./(L{3});
% 
% for i = 1:3, uIdx{i} = mIdx{i}.*(L{i} - 1); end
% [u{1}, u{2}, u{3}] = meshgrid(uIdx{:});
% 
%%


% u{1} = (mL{1} - 1).*m{1};
% u{2} = (mL{2} - 1).*m{2};
% u{3} = (mL{3} - 1).*m{3};

% shear
% k = 0.1*tand(45);
%
% u{1} = k*m{2};
% u{2} = zeros(mSize);
% u{3} = zeros(mSize);

% biaxial stretch
L(1) = 1.1;
L(2) = 1;
L(3) = 1/(L(1)*L(2));
% %
u{1} = (L(1) - 1)*m{1};
u{2} = (L(2) - 1)*m{2};
u{3} = (L(3) - 1)*m{3};

% translation with strech
% T = [0 0 0];
% L(1) = 1.1;
% L(2) = 1;
% L(3) = 1/(L(1)*L(2));
% u{1} = (L(1) - 1)*m{1} + T(1)*ones(size(m{1}));
% u{2} = (L(2) - 1)*m{2} + T(2)*ones(size(m{2}));
% u{3} = (L(3) - 1)*m{3} + T(3)*ones(size(m{3}));

uTri = triInterp(fp,u,mIdx);        % Eij on faces

vDef = triInterp(v,u,mIdx) + v;        % Eij on faces
[~, ~, cDef] = triNormals(f,vDef);

triVDef = triVolume(f,vDef);
mechanics = calculateMeanMetrics(f,v,uTri);

dm = 1;
Fij = calculateFij({u},dm,'prewitt');
FTri = triInterp(fp,Fij{1},mIdx);
ETri = calculateEijTri(FTri);

for i = 1:size(ETri,3)
    [eigVec_, eigVal_] = eig(ETri(:,:,i));
    eigVal_ = diag(eigVal_)';
    minEigVal(i,:) = [eigVec_(:,1)', eigVal_(1)];
    maxEigVal(i,:) = [eigVec_(:,3)', eigVal_(3)];   
    eigVal(i,:) = [eigVec_(:,1)', eigVal_(1), eigVec_(:,3)', eigVal_(3)];
    %             eigValW(ii,:) = [eigVec_(:,1)'*eigVal_(1,1), eigVal_(1,1), eigVec_(:,3)'*eigVal_(3,3), eigVal_(3,3)];
end

% threshold min
% min01 = eigVal(:,4) >= prctile(eigVal(:,4),50);
% max01 = eigVal(:,8) <= prctile(eigVal(:,8),50);

%         eigVal(min01,1:3) = 0;
%         eigVal(max01,5:7) = 0;
% 
        cd('C:\Users\Eyal\Dropbox\Neutrophil3D\mFiles\validation\');
dataInfo.varnames = {'x1','x2','x3','vmin1','vmin2','vmin3','Emin', 'vmax1','vmax2','vmax3','Emax'};
dataInfo.zonename = 'Cell_';
TPData = convert2Tecplot('FESurfaces', [{f},{v}],  {eigVal}, dataInfo);
mat2tecplot(TPData,'C:\Users\Eyal\Dropbox\Neutrophil3D\mFiles\validation\eigLinesValidation.plt');
%%


% eigValMag = vecMag(eigValETri,2);

figure;
view_ = {2, [90, 0], 3};
for i = 1:3
       
    subplot(1,3,i)
    h = trisurf(f,v(:,1),v(:,2),v(:,3),eigVal(:,8));
    
    
    axis(1.25*[-sphereRadius sphereRadius -sphereRadius sphereRadius -sphereRadius sphereRadius]);
    xlabel('x');
    
    ylabel('y');
    zlabel('z');

    view(view_{i})
    h.LineStyle = 'none';
    h.Parent.Color = [0 0 0];
    
    light
    axis equal
    
    if i == 2, c = colorbar('location','south'); end
    
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
    q = quiver3(fp(:,1),fp(:,2),fp(:,3),   uTri(:,1),uTri(:,2),uTri(:,3), 0);
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




    
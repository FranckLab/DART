%% runValidated Neutrophil 3D
close all; clear

sphereRadius = 32;
mSize = [256,256,256];

[f, v] = sphere_tri('ico',4,sphereRadius);
triC = [0 0 0];


v = bsxfun(@plus, v, triC);
[fn,~,fp] = triNormals(f,v);

triV = triVolume(f, v);

for i = 1:3, mIdx{i} = (1:mSize(i)) - (mSize(i)+1)/2 + triC(i); end
% [m{1}, m{2}, m{3}] = meshgrid(mIdx{:});


% uniaxial stretch
% L{3} = [linspace(1.1,1.0,mSize(3)/2), linspace(1.0,1.1,mSize(3)/2)];
%%
L{3} = [linspace(0.8,1.2,mSize(3))];

L{1} = 1;
L{2} = 1;

for i = 1:3, uIdx{i} = (mIdx{i} - (mIdx{i}-1)/2).*(L{i} - 1); end
clf; plot(uIdx{3});


%%
hold off;
[u{1}, u{2}, u{3}] = meshgrid(uIdx{:});

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
% L(1) = 1.1;
% L(2) = 0.8;
% L(3) = 1/(L(1)*L(2));
% %
% u{1} = (L(1) - 1)*m{1};
% u{2} = (L(2) - 1)*m{2};
% u{3} = (L(3) - 1)*m{3};

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

%
% mechanics = calculateMeanMetrics(f,v,uTri);

ui_nj = vecDyadic(uTri,fn);

I = eye(3,3);
F = I + (1/triV)*triIntegrate(f,v, ui_nj);

J = det(F);
C = F'*F;
U = sqrtm(C);
R = F/U;

E = 1/2*(C - I);
e = triIntegrate(f,v, ui_nj);
e = 0.5/triV*(e + e');
clc; close all;
[eigVec, eigValE] = eig(E);
eigValE = diag(eigValE);

dm = 1;
Fij = calculateFij({u},dm,'prewitt');
Eij = calculateEij(Fij);

ETri = triInterp(fp,Eij{1},dm);


for i = 1:size(ETri,3)
    eigValETri(i,:) = eig(ETri(:,:,i));
end

m = bsxfun(@minus, fp, mean(fp));
mNorm = vecNorm(m,2);
dCdt = [0 0 1];
dCdtNorm = dCdt/norm(dCdt);

%   binTri = sum(bsxfun(@times, m , dCdtNorm),2);
%   bins{i}(:,j) = linspace(min(binTri),max(binTri), nBins)
binTri = sum(bsxfun(@times, mNorm, dCdtNorm),2);

nBins = 50;
bins = linspace(-1,1, nBins + 2);
bins = bins(2:end-1);
[~,binInd] = histc(binTri, bins);



for i = 1:nBins
    eigValETri_ = eigValETri(binInd == i,:);
    
    if isempty(eigValETri_),
        maxEig(i) = nan;
        minEig(i) = nan;
    else
        maxEig(i) = max(eigValETri_(:,3));
        minEig(i) = min(eigValETri_(:,1));
    end
    
end


%%


% eigValMag = vecMag(eigValETri,2);

figure;
view_ = {2, [90, 0], 3};
for i = 1:3
       
    subplot(1,3,i)
    h{1} = trisurf(f,v(:,1),v(:,2),v(:,3),vecMag(eigValETri,2));
   
    
    
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




    
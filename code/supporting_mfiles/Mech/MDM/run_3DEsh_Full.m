%% runValidated Neutrophil 3D
close all; clear

load('C:\Users\Eyal\Dropbox\Neutrophil3D\mFiles\validateEshelby\resultsEsh_xpos_uneqyzneg.mat');
u0 = u;

 savepath = 'C:\Users\Eyal\Google Drive\Research\Publications\14 Neutrophil3D\manuscript\figures\SI\eshelby\141208_eshelby';
%%
u = u0{1};

sphereRadius = 32;
mSize = 129*[1 1 1];

[f, v] = sphere_tri('ico',4,sphereRadius);
triV = triVolume(f, v);
[fn, ~, fp] = triNormals(f,v);

for i = 1:3, mIdx{i} = (1:mSize(i)) - (mSize(i)+1)/2; end
[m{1}, m{2}, m{3}] = meshgrid(mIdx{:});
BW = sqrt(m{1}.^2 + m{2}.^2 + m{3}.^2) <= 32;


minv = floor(min(v));
maxv = ceil(max(v));


u{4} = sqrt(u{1}.^2 + u{2}.^2 + u{3}.^2);
u{5} = bwdist(BW);

minv = floor(min(v)*3);
maxv = ceil(max(v)*3);

nPoints = 30*1000;
xS = bsxfun(@times, rand(nPoints,3), maxv - minv);
xS = bsxfun(@plus, xS, minv + 64);

%         for k = 1:5
%             u{k} = permute(u{k},[2 1 3]);
%         end
uS           = triInterp(xS,u,[1 1 1]);
uS(:,4) = sqrt(uS(:,1).^2 + uS(:,2).^2 + uS(:,3).^2);
xS = xS - 64;
%%
        distanceThreshold = 100;

        idx01 = uS(:,5) ~= 0;
      

X = xS(idx01,1);
Y = xS(idx01,2);
Z = xS(idx01,3);
U = uS(idx01,1);
V = uS(idx01,2);
W = uS(idx01,3);
C = uS(idx01,4);
%         C = colormap(jet(length(uS(:,4))));

%         close all
figure;
hold all;
hc = coneplot(X,Y,Z,U,V,W,0.03,'nointerp');
fvc = repmat(C(:)',[42 1]);
set(hc, 'facecolor', 'flat', 'facevertexcdata', fvc(:))
hc.EdgeColor = 'none';
    hc.AmbientStrength = 0.6;
        hc.DiffuseStrength = 0.75;
        hc.SpecularStrength = 0.4;
        hc.SpecularExponent = 3;
colormap(hot)
freezeColors
hs = trisurf(f,v(:,1),v(:,2),v(:,3),'facecolor',[85 85 85]/255,'edgecolor','none');
  hs.AmbientStrength = 0.25;
        hs.DiffuseStrength = 0.50;
        hs.SpecularStrength = 0.4;
        hs.SpecularExponent = 3;
axis image;
hl = light;
lightangle(hl,160, 20)

lighting gouraud

h = gcf;
set(gca,'color','none','Visible','off')
set(h,'color','none');
set(h, 'InvertHardCopy', 'off')
set(h, 'PaperUnits','inches')
set(h, 'PaperSize',[3 4.5])
set(h, 'Paperposition', [3 3 3 4.5]);
set(h, 'Units','inches','Position', [3, 3, 3 4.5]);
%%
view([132.61 6.5]);
savename = fullfile(savepath,'u_3D_cone.png');
export_fig(h, savename,  '-png', '-opengl', '-r1000');

view([180 0]);
savename = fullfile(savepath,'u_3D_cone_Front.png');
export_fig(h, savename,  '-png', '-opengl', '-r1000');

view([90 0]);
savename = fullfile(savepath,'u_3D_cone_Side.png');
export_fig(h, savename,  '-png', '-opengl', '-r1000');
%%
[~, ~, fp] = triNormals(f,v);
   uTri           = triInterp(fp + 129/2,u,[1 1 1]);
% uTri = uS(:,1:3);
% idx01 = uTri(:,5) ~= 0;
idx01 = 1:size(uTri,1);
        uTriMag = vecMag(uTri,2);
        
%         idx = 1:length(fp);
%         idx01 = round(linspace(1,length(fp),5000));

        X = fp(idx01,1);
        Y = fp(idx01,2);
        Z = fp(idx01,3);
        
        U = uTri(idx01,1);
        V = uTri(idx01,2);
        W = uTri(idx01,3);
        C = uTriMag(idx01);
        
       
        close all
        figure;
        hold all;
        hc = coneplot(X,Y,Z,U,V,W,0.03,'nointerp');
        fvc = repmat(C(:)',[42 1]);
        set(hc, 'facecolor', 'flat', 'facevertexcdata', fvc(:))
        hc.EdgeColor = 'none';
        hc.AmbientStrength = 0.25;
        hc.DiffuseStrength = 0.50;
        hc.SpecularStrength = 0.4;
        hc.SpecularExponent = 3;
        colormap(hot)
        caxis([0 4.25])
        
        freezeColors
        hs = trisurf(f,v(:,1),v(:,2),v(:,3),'facecolor',[42 42 42]/255,'edgecolor','none');
        hs.AmbientStrength = 0.25;
        hs.DiffuseStrength = 0.50;
        hs.SpecularStrength = 0.4;
        hs.SpecularExponent = 3;
        axis image;
        hl = light;
        lightangle(hl,160, 20)
        view([132.61 6.5]);
        
        lighting gouraud
        
        h = gcf;
        set(gca,'color','none','Visible','off')
        set(h,'color','none');
        set(h, 'InvertHardCopy', 'off')
        set(h, 'PaperUnits','inches')
        set(h, 'PaperSize',[3 4.5])
        set(h, 'Paperposition', [3 3 3 4.5]);
        set(h, 'Units','inches','Position', [3, 3, 3 4.5]);
        
        %
       
        savename = fullfile(savepath,'u_3DConeSurface.png');
        export_fig(h, savename,  '-png', '-opengl', '-r1000');
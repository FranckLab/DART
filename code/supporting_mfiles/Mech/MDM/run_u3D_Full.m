
% Folder Directory

clc; clear;% close all;
% path = 'F:\Neutrophil3D\Data\';
% path = '\\THEYETI\Neutrophil3D\Neutrophil3D\Data\';
% path = 'E:\Neutrophil3D\Data';
path = 'C:\Users\Eyal\Google Drive\Research\Publications\14 Neutrophil3D\data';
folderDir{01} = 'chemotaxis\Naive\neutrCollagen_1uM';
% folderDir{02} = 'chemotaxis\Naive\neutrCollagen_1uM_cell3';
% folderDir{03} = 'chemotaxis\Naive\neutrCollagen_100nM';

%%

for i = 1:length(folderDir)
    cd(fullfile(path,folderDir{i}));
    disp([num2str(i),'  ', folderDir{i}]);
    
    
    load('resultsDVC.mat','u','um2vxl','dm');
    load('resultsCell.mat');
    
    nTime   = min([length(u), length(f) - 1]);       % number of DVC increments points
    % --------------------------------------------------------n---------------------------------------
    % -----------------------------------------------------------------------------------------------
    
    
    u0 = u;
    u         = changeUnits(u0,um2vxl);
    uTri      = cell(nTime,1);
    
    for j = 5%%:nTime
        
        %% calculate 3D tri data
        
        [u_, uMean(i,:), uSTD(i,:), uPValue(i,:)] = removeNoise_Neutrophil3D(u{j},2);
        %         u_ = u{j};
        u{j} = u_;
        [fn, ~, fp] = triNormals(f{j},v{j});
        
        
        for k = 1:length(u{j})
            u{j}{k} = u{j}{k} - mean(u{j}{k}(:));
        end
        
        
        
        uTri{j}(:,4)      = vecMag(uTri{j},2);
        
        tri = triangulation(f{j},v{j}(:,1),v{j}(:,2),v{j}(:,3));
        
        v_ = bsxfun(@rdivide, v{j}, 4*um2vxl);
        ch = convhulln(v_);
        
        minv = floor(min(v_));
        maxv = ceil(max(v_));
        
        
        [m{1},m{2},m{3}] = ndgrid(minv(2):1:maxv(2), minv(1):1:maxv(1), minv(3):1:maxv(3));
        m{4} = [m{2}(:) m{1}(:) m{3}(:)];
        in = intriangulation(v_,ch,m{4});
        
        uSize = size(u{1}{1});
        BW = zeros(size(u{1}{1}));
        BW(minv(2):maxv(2), minv(1):maxv(1), minv(3):maxv(3)) = reshape(in,size(m{1}));
        % %         BW{j,1} = bwdist(BW_);
        %         BW{j,1} = bwdist(BW_);
        
        u{j}{4} = sqrt(u{j}{1}.^2 + u{j}{2}.^2 + u{j}{3}.^2);
        u{j}{5} = bwdist(BW);
        
        meanV = mean(v{j});
        v_ = bsxfun(@minus, v{j}, meanV);
        
        minv = floor(min(v_)*3);
        minv(3) = -12;
        maxv = ceil(max(v_)*3);
        
        nPoints = 30*1000;
        xS{j} = bsxfun(@times, rand(nPoints,3), maxv - minv);
        xS{j} = bsxfun(@plus, xS{j}, minv + meanV);
        
        %         for k = 1:5
        %             u{j}{k} = permute(u{j}{k},[2 1 3]);
        %         end
        uS{j}           = triInterp(xS{j},u{j}(:,1:5),dm);
        uS{j}(:,4) = sqrt(uS{j}(:,1).^2 + uS{j}(:,2).^2 + uS{j}(:,3).^2);
        %%
        distanceThreshold = 20;
        
        idx01 = uS{j}(:,5) <= distanceThreshold;
        
        X = xS{j}(idx01,1);
        Y = xS{j}(idx01,2);
        Z = xS{j}(idx01,3);
        U = uS{j}(idx01,1);
        V = uS{j}(idx01,2);
        W = uS{j}(idx01,3);
        C = uS{j}(idx01,4);
        %         C = colormap(jet(length(uS{j}(:,4))));
        %%
        close all
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
        caxis([0 1])
        freezeColors
        
        
        hs = trisurf(f{j},v{j}(:,1),v{j}(:,2),v{j}(:,3),'facecolor',[60 60 60]/255,'edgecolor','none');
        hs.AmbientStrength = 0.40;
        hs.DiffuseStrength = 0.750;
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
         savepath = 'C:\Users\Eyal\desktop';
%         savepath = 'C:\Users\Eyal\Google Drive\Research\Publications\14 Neutrophil3D\manuscript\figures\mapping\141205';
        savename = fullfile(savepath,'u3DCone1.png');
        export_fig(h, savename,  '-png', '-opengl', '-r1000');
        %         print('-dpng','-opengl','-r1000',savename)
        %%
        
        
        uTri{j}           = triInterp(fp,u_,dm);
        uTriMag = vecMag(uTri{j},2);
        
%         idx = 1:4:length(fp);
        idx = round(linspace(1,length(fp),7500));

        X = fp(idx,1);
        Y = fp(idx,2);
        Z = fp(idx,3);
        
        U = uTri{j}(idx,1);
        V = uTri{j}(idx,2);
        W = uTri{j}(idx,3);
        C = uTriMag(idx);
        
        
        
        
        
        close all
        figure;
        hold all;
        hc = coneplot(X,Y,Z,U,V,W,0.05,'nointerp');
        fvc = repmat(C(:)',[42 1]);
        set(hc, 'facecolor', 'flat', 'facevertexcdata', fvc(:))
        hc.EdgeColor = 'none';
        hc.AmbientStrength = 0.25;
        hc.DiffuseStrength = 0.50;
        hc.SpecularStrength = 0.4;
        hc.SpecularExponent = 3;
        colormap(hot)
        caxis([0 1])
        
        freezeColors
        hs = trisurf(f{j},v{j}(:,1),v{j}(:,2),v{j}(:,3),'facecolor',[42 42 42]/255,'edgecolor','none');
        hs.AmbientStrength = 0.6;
        hs.DiffuseStrength = 0.75;
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
        
         savepath = 'C:\Users\Eyal\desktop';
%         savepath = 'C:\Users\Eyal\Google Drive\Research\Publications\14 Neutrophil3D\manuscript\figures\mapping\141205';
        savename = fullfile(savepath,'u3DConeSurface.png');
        export_fig(h, savename,  '-png', '-opengl', '-r1000');
        %         print('-dpng','-opengl','-r1000',savename)
        %%
        
    
    end
    
    dataInfo.varnames = {'x1','x2','x3','ux','uy','uz','uMag','BW'};
    dataInfo.zonename = 'u3D_';
    TPData = convert2Tecplot('lines', xS, uS, dataInfo);
    mat2tecplot(TPData,'u3DFull.plt');
    
    %     dataInfo.varnames = {'x1','x2','x3','ux','uy','uz','uMag','BW'};
    %     dataInfo.zonename = 'u3D_';
    %     TPData = convert2Tecplot('cubes', um2vxl*4, u, dataInfo);
    %     mat2tecplot(TPData,'u3DFull.plt');
    
    %     dataInfo.varnames = {'x1','x2','x3',''};
    %     dataInfo.zonename = 'Cell_';
    %     TPData = convert2Tecplot('FESurfaces', [f(1:end-1),v(1:end-1)],  [mEmTri, dudtTri], dataInfo);
    %     mat2tecplot(TPData,'cell.plt');
    
end

% Folder Directory

clc; clear; %close all;
% path = 'F:\Neutrophil3D\Data\';
% path = '\\THEYETI\Neutrophil3D\Neutrophil3D\Data\';
% path = 'E:\Neutrophil3D\Data';
path = 'C:\Users\Eyal\Dropbox\Neutrophil3D\data';


folderDir{01} = 'chemotaxis\LPS\neutr_10_20_13_activ_cell2';
folderDir{02} = 'chemotaxis\LPS\neutr_10_20_13_activ_cell3';
folderDir{03} = 'chemotaxis\LPS\neutr_10_20_13_activ_cell4';
folderDir{04} = 'chemotaxis\Naive\neutrCollagen_1uM';
folderDir{05} = 'chemotaxis\Naive\neutrCollagen_1uM_cell3';
folderDir{06} = 'chemotaxis\Naive\neutrCollagen_100nM';
folderDir{07} = 'chemokinesis\LPS\neutr_kineses_Apr7cell2';
folderDir{08} = 'chemokinesis\LPS\neutr_kineses_Apr23cell1';
folderDir{09} = 'chemokinesis\LPS\neutr_kineses_Apr23cell2';
folderDir{10} = 'chemokinesis\LPS\neutr_kineses_May8cell1';
folderDir{11} = 'chemokinesis\LPS\neutr_kineses_May8cell3';
folderDir{12} = 'chemokinesis\Naive\neutr_kineses_Apr7cell1';
folderDir{13} = 'chemokinesis\Naive\neutr_kineses_Apr11cell1';
folderDir{14} = 'chemokinesis\Naive\neutr_kineses_Apr22cell2';
folderDir{15} = 'chemokinesis\Naive\neutr_kineses_May8cell2';

%%

for i = 4%:length(folderDir)
    cd(fullfile(path,folderDir{i}));
    disp([num2str(i),'  ', folderDir{i}]);
    
    load('resultsDVC.mat','u','um2vxl','dm');
    load('resultsCell.mat');
    u0 = u;
    u         = changeUnits(u0,um2vxl);
    
    nTime   = min([length(u), length(f) - 1]);       % number of DVC increments points
    
    % -----------------------------------------------------------------------------------------------
    % -----------------------------------------------------------------------------------------------
    
    
    eigVal = cell(nTime,1);
    
    for j = 5%:nTime
        
        %         fp = triCenters(f{j},v{j});
        [fn, ~, fp] = triNormals(f{j},v{j});
        Fij = calculateFij(u(j),4,'optimal9');
        FTri = triInterp(fp,Fij{1},dm);
        ETri = calculateEijTri(FTri);
        
        
        uTri = triInterp(fp,u{5}, dm);
        
        
        
        eigVal{j,1} = zeros(size(ETri,3),8);
        eigValNorm = zeros(size(ETri,3),3);
        L = zeros(size(FTri,3),3);
        
        for ii = 1:size(ETri,3)
            [eigVec_, eigVal_] = eig(ETri(:,:,ii));
            eigVal_ = diag(eigVal_)';
            eigVal{j,1}(ii,:) = [eigVec_(:,1)', eigVal_(1), eigVec_(:,3)', eigVal_(3)];
            
            eigValNorm(ii,:) = eigVal_/norm(eigVal_);
            
            %             L1 >= L2 >= L3
            F_ = FTri(:,:,ii) - eye(3);
            L_ = eig(F_);
            [~, sortIdx] = sort(real(L_),'descend');
            L(ii,:) = L_(sortIdx);
            
            n = fn(ii,:);
            R =  vrrotvec([0 0 1],n);
            
        end
        
        %%
        
        
        
        
        %%
        f_ = f{5};
        v_ = v{5};
        t = triangulation(f_,v_);
        f0 = 15148;
        
        f1 = neighbors(t,f0)';
        f2 = cell2mat(vertexAttachments(t,f_(f0,:)')')';
        f2(f2 == f1(1)) = [];
        f2(f2 == f1(2)) = [];
        f2(f2 == f1(3)) = [];
        f2(f2 == f0) = [];
        
        f3 = [f0; f1; f2];
        L(f3,:)
        
        
        
        FTri(:,:,f3)
        mean(ETri(:,:,f3),3)
        
        
        figure;
        
        hold on;
        trisurf(t,'facecolor',[85 85 85]/255,'edgecolor','none');
        trimesh(f_(f0,:),v_(:,1),v_(:,2),v_(:,3),'facecolor',[255 0 0]/255);
        trimesh(f_(f1,:),v_(:,1),v_(:,2),v_(:,3),'facecolor',[0 255 0]/255);
        trimesh(f_(f2,:),v_(:,1),v_(:,2),v_(:,3),'facecolor',[0 0 255]/255);
        %         trimesh(f_(idx,:),v_(:,1),v_(:,2),v_(:,3),'facecolor',[255 0 0]/255);
        
        
        axis image;
        hl = light;
        lightangle(hl,160, 20)
        lighting gouraud
        view([132.61 6.5]);
        hold off
        
        %%
        
        
        node01      =  find(~any(imag(L),2));
        nodeRep01   =   node01(sum(L(node01,:) > 0,2) == 3);
        nodeAtt01   =   node01(sum(L(node01,:) > 0,2) == 0);
        saddleRep01 =   node01(sum(L(node01,:) > 0,2) == 2);
        saddleAtt01 =   node01(sum(L(node01,:) > 0,2) == 1);
        
        
        spiral01 = find(any(imag(L),2));
        spiralRep01 = spiral01(sum( real(L(spiral01,:)) > 0,2) == 3);
        spiralAtt01 = spiral01(sum( real(L(spiral01,:)) > 0,2) == 0);
        spiralSaddleRep01 = spiral01(sum( real(L(spiral01,:)) > 0,2) == 2);
        sprialSaddleAtt01 = spiral01(sum( real(L(spiral01,:)) > 0,2) == 1);
        
        
        type = zeros(size(L(:,1)));
        type(nodeRep01) = 1;
        type(nodeAtt01) = 2;
        type(saddleRep01) = 3;
        type(saddleAtt01) = 4;
        
        type(spiralRep01) = 5;
        type(spiralAtt01) = 6;
        type(spiralSaddleRep01) = 7;
        type(sprialSaddleAtt01) = 8;
        f_ = f{5};
        v_ = v{5};
        
        figure;
        hold on;
        trisurf(f_,v_(:,1),v_(:,2),v_(:,3),type,'edgecolor','none');
        %         trimesh(f_(f1,:),v_(:,1),v_(:,2),v_(:,3),'facecolor',[0 255 0]/255);
        %         trimesh(f_(f2,:),v_(:,1),v_(:,2),v_(:,3),'facecolor',[0 0 255]/255);
        %         trimesh(f_(idx,:),v_(:,1),v_(:,2),v_(:,3),'facecolor',[255 0 0]/255);
        cm(1,:) = [96 57 19]/255;
        cm(2,:) = [255 102 51]/255;
        cm(3,:) = [0 255 0]/255;
        cm(4,:) = [255 51 153]/255;
        cm(5,:) = [51 53 204]/255;
        cm(6,:) = [0 0 255]/255;
        cm(7,:) = [251 251 55]/255;
        cm(8,:) = [255 0 0]/255;
        colormap(cm)
        
        axis image;
        hl = light;
        lightangle(hl,160, 20)
        lighting flat
        view([132.61 6.5]);
        hold off
        
        set(gca,'color',[0 0 0])
        
        
        
        figure;
        hold on;
        trisurf(f_,v_(:,1),v_(:,2),v_(:,3),vecMag(uTri,2),'edgecolor','none');
        axis image;
        hl = light;
        lightangle(hl,160, 20)
        lighting flat
        view([132.61 6.5]);
        hold off
        colormap(hot);
        set(gca,'color',[0 0 0])
        %%
        %         real
        x = zeros(size(L,1),1);
        y = zeros(size(L,1),1);
        z = zeros(size(L,1),1);
        
        
        
        % equation 1: for all real eigenvales
        allReal = ~any(imag(L),2);
        x(allReal) = L(allReal,2);
        y(allReal) = min([L(allReal,1) - L(allReal,2), L(allReal,2) - L(allReal,3)],[],2);
        z(allReal) = 2*L(allReal,2) - (L(allReal,1) + L(allReal,3));
        
        % equation 2: for all real eigenvales
        imag01 = imag(L) ~= 0;
        c = zeros(size(L,1),1);
        a = zeros(size(L,1),1);
        b = zeros(size(L,1),1);
        
        c01 = ((imag(L) == 0) & repmat(~all(imag(L) == 0,2),1,3));
        a01 = ((imag(L) ~= 0) & repmat(~all(imag(L) == 0,2),1,3));
        b01 = ((imag(L) ~= 0) & repmat(~all(imag(L) == 0,2),1,3));
        
        c(any(c01,2)) = L(c01);
        a(any(a01,2)) = real(L(a01 & imag(L) > 0));
        b(any(b01,2)) = imag(L(b01 & imag(L) < 0));
        
        x(~allReal) = a(~allReal);
        y(~allReal) = b(~allReal);
        z(~allReal) = a(~allReal) - c(~allReal);
        
        r = sqrt(x.^2 + y.^2 + z.^2);
        x = x./r;
        y = y./r;
        z = z./r;
        
        u_minus_x = x <= 0;
        alpha = asin(y);
        beta = atan2(z,x) + pi/2*u_minus_x.*sign(z);
        
    end
end
%%


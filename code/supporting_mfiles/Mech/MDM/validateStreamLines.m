
% Folder Directory

clc; clear; %close all;
% path = 'F:\Neutrophil3D\Data\';
% path = '\\THEYETI\Neutrophil3D\Neutrophil3D\Data\';
% path = 'E:\Neutrophil3D\Data';
path = 'C:\Users\Eyal\Dropbox\Neutrophil3D\data';


% folderDir{01} = 'chemotaxis\LPS\neutr_10_20_13_activ_cell2';
% folderDir{02} = 'chemotaxis\LPS\neutr_10_20_13_activ_cell3';
% folderDir{03} = 'chemotaxis\LPS\neutr_10_20_13_activ_cell4';
% folderDir{04} = 'chemotaxis\Naive\neutrCollagen_1uM';
% folderDir{05} = 'chemotaxis\Naive\neutrCollagen_1uM_cell3';
folderDir{01} = 'chemotaxis\Naive\neutrCollagen_100nM';
% folderDir{07} = 'chemokinesis\LPS\neutr_kineses_Apr7cell2';
% folderDir{08} = 'chemokinesis\LPS\neutr_kineses_Apr23cell1';
% folderDir{09} = 'chemokinesis\LPS\neutr_kineses_Apr23cell2';
% folderDir{10} = 'chemokinesis\LPS\neutr_kineses_May8cell1';
% folderDir{11} = 'chemokinesis\LPS\neutr_kineses_May8cell3';
% folderDir{12} = 'chemokinesis\Naive\neutr_kineses_Apr7cell1';
% folderDir{13} = 'chemokinesis\Naive\neutr_kineses_Apr11cell1';
% folderDir{14} = 'chemokinesis\Naive\neutr_kineses_Apr22cell2';
% folderDir{15} = 'chemokinesis\Naive\neutr_kineses_May8cell2';

%%

n = 512;
[mBW{1},mBW{2}] = meshgrid(-(n - 1)/2:(n - 1)/2 ,-(n - 1)/2:(n - 1)/2);
BW = sqrt(mBW{1}.^2 + mBW{2}.^2) < (n - 1)/2;
nanBW = nan(n,n);
nanBW(BW) = 1;


for i = 1
    cd(fullfile(path,folderDir{i}));
    disp([num2str(i),'  ', folderDir{i}]);
    
    load('resultsDVC.mat','u','um2vxl','dm');
    load('resultsCell.mat');
    
    nTime   = min([length(u), length(f) - 1]);       % number of DVC increments points
    
    % -----------------------------------------------------------------------------------------------
    % -----------------------------------------------------------------------------------------------
    
    u0 = u;
    u         = changeUnits(u0,um2vxl);
    uTri      = cell(nTime,1);
    TP = cell(nTime,1);
    TPMap = cell(nTime,1);
    
    for j = 1%:nTime
        
        %% calculate 3D tri data
        [fn, ~, fp] = triNormals(f{j},v{j});
        R = mean(vecMag(bsxfun(@minus, fp, mean(fp)),2));        

        uTri{j}           = triInterp(fp,u{j},dm);
        uTri{j}(:,4)      = vecMag(uTri{j},2);
        
%         uNan = zeros(size(uTri{j}));
%         uNan([2300,3000,4000,2500,1000],:) = uTri{j}([2300,3000,4000,2500,1000],:);
%          uTri{j} = uNan;
        
        uN = bsxfun(@times, fn,sum( uTri{j}(:,1:3).*fn ,2));
        uT = uTri{j}(:,1:3) - uN;
        
        uN = sum( uTri{j}(:,1:3).*fn,2);
        %         uP = vecMag(uP,2);
        
        [Fij,~, gradU] = calculateFij(u0(j),4,'optimal9');
        FTri = triInterp(fp,Fij{1},dm);
        I = repmat(eye(3,3), 1, 1, size(FTri,3));
        gradUTri = FTri - I;
        gradUTri = sqrt(sum(sum(gradUTri.*gradUTri,1),2));
        gradUTri = gradUTri(:);
        
        %         {'x1','x2','x3','u1','u2','u3','|u|','uP','uN','gradU'};
        TP{j} = [uTri{j}, uT, vecMag(uT,2), uN, gradUTri];
        
        %% mapping
        
        m = bsxfun(@minus, fp, mean(fp));
        % [phi, theta, r]
        [long,latm,r] = cart2sph(m(:,1),m(:,2),m(:,3));
        
        figure;
        axesm ('mollweid', 'Frame', 'on', 'Grid', 'on');
        mstruct = gcm;
        z = r;
        [x,y,z,savepts] = mfwdtran(mstruct,latm*180/pi,long*180/pi,z,'linem');
        
        nanIdx = isnan(x);
         x(isnan(x)) = [];
         y(isnan(y)) = [];
         z(isnan(z)) = [];
         
         
         clear uTSph
         lat = pi/2 - latm;
         phi = long;
         theta = lat;
        uTSph(:,1) = sin(theta).*cos(phi).*uT(:,1) + sin(theta).*sin(phi).*uT(:,2) + cos(theta).*uT(:,3);
        uTSph(:,2) = -(cos(theta).*cos(phi).*uT(:,1) + cos(theta).*sin(phi).*uT(:,2) - sin(theta).*uT(:,3));
        uTSph(:,3) = -sin(phi).*uT(:,1) + cos(phi).*uT(:,2);
%         quiverm(latm*180/pi,long*180/pi,uTSph(:,2),uTSph(:,3))
  
%         theta = lat;
%         conv = 0;
%         for ii = 1:1000
%             theta0 = theta;
%             theta = theta - (2*theta + sin(2*theta) - pi*sin(lat))./(2 + 2*cos(2*theta));
%         end
%         
%         rotation = 0;
%         
%         longRot = long;
%         longRot( long < 0) = longRot( long < 0) + 2*pi;
%          longRot = longRot - pi;        
%         x = 1*2*sqrt(2)/pi.*(longRot).*cos(theta);
%         y =  1*sqrt(2).*sin(theta);

        [mx, my] = meshgrid(linspace(min(x),max(x),n), linspace(min(y),max(y),n));
            
        F = scatteredInterpolant([x,y], uN);
        uNMap = F(mx,my);
        
        uTMap = cell(1,4);
        for ii = 1:3
            F = scatteredInterpolant([x,y], uTSph(:,ii));
            uTMap{ii} = F(mx,my);
        end
        uTMap{4} =  sqrt(uTMap{2}.^2 + uTMap{3}.^2);
        
        %         {'x1','x2','x3','uT1','uT2','uT3', 'uT4','uNMap'}
        TPMap{j} = [uTMap(1), uTMap(2), uTMap(3), uTMap(4), {uNMap}, {BW}];
        
        
    end
    
    
    dataInfo.varnames = {'x1','x2','x3','u1','u2','u3','|u|','uT1','uT2','uT3','|uT|','uN','gradU'};
    dataInfo.zonename = 'Cell_';
    TPData = convert2Tecplot('FESurfaces', [f(1:end-1),v(1:end-1)],  TP(1), dataInfo);
    mat2tecplot(TPData,'3DTri.plt');
    
    dataInfo.varnames = {'x1','x2','x3','uT1','uT2','uT3', 'uT4','uNMap','BW'};
    dataInfo.zonename = 'Cell_';
    TPData = convert2Tecplot('surfaces', {mx, my},  TPMap(1), dataInfo);
    mat2tecplot(TPData,'2DMap.plt');
    
    %         dataInfo.varnames = {'x1','x2','x3','u1','u2','u3','uMag' ,'uR','uN','uT'};
    %         dataInfo.zonename = 'Cell_';
    %         TPData = convert2Tecplot('FESurfaces', [f(1:end-1),v(1:end-1)],  uTri, dataInfo);
    %         mat2tecplot(TPData,'u.plt');
    %
    %
    %     dataInfo.varnames = {'x1','x2','x3','u1','u2','u3', 'uR','uTheta','uPhi','uMag', 'BW'};
    %     dataInfo.zonename = 'Cell_';
    %     TPData = convert2Tecplot('surfaces', {mx, my},  uMapSph, dataInfo);
    %     mat2tecplot(TPData,'uMap.plt');
    %
    %     dataInfo.varnames = {'x1','x2','x3','graduMag'};
    %     dataInfo.zonename = 'Cell_';
    %     TPData = convert2Tecplot('FESurfaces', [f(1:end-1),v(1:end-1)],  gradUTri, dataInfo);
    %     mat2tecplot(TPData,'gradU.plt');
    %
    %
    %     dataInfo.varnames = {'x1','x2','x3','gradUMagMap','BW'};
    %     dataInfo.zonename = 'Cell_';
    %     TPData = convert2Tecplot('surfaces', {mx, my},  gradUMap, dataInfo);
    %     mat2tecplot(TPData,'gradUMap.plt');
    %
    
    
end

%%
% figure
% 
% [~,h] = contourf(mx,my,sqrt(uTMap{2}.^2 + uTMap{3}.^2).*nanBW,20);  colorbar; colormap parula
% set(h,'linestyle','none');
% 
% hold on;
% [mx, my] = meshgrid(linspace(min(x),max(x),n), linspace(min(y),max(y),n)); 
% streamslice(mx,my,uTMap{3}.*nanBW,uTMap{2}.*nanBW,8)
% 
% dq = 12;
% qmx = mx(:);
% qmy = my(:);
% qmx = qmx(1:dq:end);
% qmy = qmy(1:dq:end);
% quTMap2 = uTMap{2}.*nanBW;
% quTMap3 = uTMap{3}.*nanBW;
% quTMap2 = quTMap2(:);
% quTMap3 = quTMap3(:);
% quTMap2 = quTMap2(1:dq:end);
% quTMap3 = quTMap3(1:dq:end);
% quiver(qmx,qmy,quTMap3,quTMap2,1);
% 
% hold off;
% axis image
% 
% 
% figure;
% dq = 1;
% h = trisurf(f{j},v{j}(:,1),v{j}(:,2),v{j}(:,3),vecMag(uTSph,2));
% set(h,'linestyle','none');
% qfp = fp(1:dq:end,:);
% quT = uT(1:dq:end,:);
% 
% hold on;
% q =quiver3(qfp(:,1),qfp(:,2),qfp(:,3),quT(:,1),quT(:,2),quT(:,3),3,'color','k');
% hold off
%                 a = gca;
%                 a.CLim = [-.2 .2];
% axis image;
% xlabel('x'); ylabel('y'); zlabel('z');
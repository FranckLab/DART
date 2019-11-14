
% Folder Directory

clc; clear; close all;
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
mechanics = cell(length(folderDir),1);

for i = 1:length(folderDir)
    cd(fullfile(path,folderDir{i}));
    disp([num2str(i),'  ', folderDir{i}]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('resultsCell.mat');
    nTime   = length(f);       % number of DVC increments points
    
    CA = zeros(nTime,3);
    CV = zeros(nTime,3);
    for j = 1:nTime
        [CA(j,:), CV(j,:)] = triCentroid(f{j},v{j});
    end
    
    dC = diff(CA);
    
    
    
    ma = msdanalyzer(2, 'µm', 's')
    
    
    %%
    
    
    
%     
%     N = nTime;
%     
%     
%     for k = 1:N
%         for j = 1:N-k
%             d = CA(k + j,:) - CA(k,:)
%             xik(j,k,:) = d.*d;
%         end
%     end
%     
%     % non overlapping
%     xiBar = zeros(N-1,3);
%     for j = 1:N
%         ni = floor((N - 1)/j)
%         for k = 0:(ni - 1)
%             xiBar_ = 1/ni*xik(j,1+j*k,:);
%            xiBar(j,:) =  xiBar(j,:) + xiBar_(:)';
%         
%         end
%     end
%     
%     % non overlapping
%     xiBar = zeros(N-1,3);
%     for j = 1:N
%         ni = 1
%         for k = 0:(N - ni)
%             xiBar_ = 1/ni*xik(j,1+j*k,:);
%            xiBar(j,:) =  xiBar(j,:) + xiBar_(:)';
%         
%         end
%     end
%     
%     
%     
%     % overlapping
%     ni2 = squeeze(sum(xik > 0))
%     xiBar2 = squeeze(sum(xik,2))./ni2
%         
%     for ii = 1:N
%     
%         
%     ni = floor((N - 1)/ii);
%     
% %     1/n*
%     
%     
%     
%     end
%     
%     
%     totalDistance = (vecMag(dC,2));
%     
%     phi = acos(vecCos(dC(1:end-1,:),dC(2:end,:),2));
%     
%     totalTime = 2*(numel(phi) - 1);
%     t0 = 0:2:totalTime;
%     ii = 1;
%     
%     clear phiRMS
%     dt = 1e-5;
%     dt_ = 0:dt:dt*1000;
%     for idt = dt_
%         
%         t = 0:idt:dt*1000;
%         phi_ = interp1(t0, phi', t,'spline');
%         phiRMS(ii) = sqrt(1/numel(phi_)*sum(phi_.^2));
%         ii = ii + 1;
%     end
% %     
% %     figure; plot(t0,phi);
% %     hold on;
% %     plot(t,phi_,'r:');
% %     
%     
%     phiRMS2 = phiRMS.^2;
%     
% %     persistenceTime = mean(dt_(2:100)./phiRMS2(2:100))
%     figure; plot( dt_,1./phiRMS2','.');
%     
%     
%     figure; plot(dt_(2:100),2*dt_(2:100)./phiRMS2(2:100))
    
%     %%
%     xDistance = (abs(dC(:,1)));
%     histVar{i} = xDistance./totalDistance; 
    
    
    
    
    
    
    
    dt = 2;
    centroid.CA = CA;
    centroid.CV = CV;
    centroid.dCAdt = diff(CA)/dt;
    centroid.dCVdt = diff(CV)/dt;
    
    
    
    
    
    save('resultsCell.mat','f','v','centroid');
    
    
    
end





%%
for i = 1:length(folderDir)
    cd(fullfile(path,folderDir{i}));
    disp([num2str(i),'  ', folderDir{i}]);
    

    load('resultsCell.mat', 'centroid');
        centroid_{i,1} = centroid;
        
    
    
end

clear centroid;
centroid = centroid_;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsCentroid.mat','centroid');




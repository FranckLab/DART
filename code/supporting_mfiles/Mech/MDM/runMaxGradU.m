
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

maxGradU = cell(length(folderDir),1);
for i = 1:length(folderDir)
    cd(fullfile(path,folderDir{i}));
    disp([num2str(i),'  ', folderDir{i}]);
    
    load('resultsDVC.mat','u','um2vxl','dm');
    load('resultsCell.mat');
    
    nTime   = min([length(u), length(f) - 1]);       % number of DVC increments points
    
    % -----------------------------------------------------------------------------------------------
    % -----------------------------------------------------------------------------------------------
    
    
    gradU = cell(nTime,1);
    
    for j = 1:nTime
        
        fp = triCenters(f{j},v{j});
        
         u_ = removeNoise_Neutrophil3D(u{j},2);
         
        [Fij,~, gradU] = calculateFij({u_},4,'optimal9');
        
        FTri = triInterp(fp,Fij{1},dm);
        
        nTri = size(FTri,3);
        I = repmat(eye(3,3), 1, 1, nTri);
        
        gradUTri = FTri - I;
        gradUTri = reshape(sqrt(sum(sum(gradUTri.*gradUTri,1),2)),[],1)';
        
        maxGradUTri{i}(j,1) = max(gradUTri);
        maxGradU{i}(j,1) = max(gradU{1}(:));
    end
  
    
end

save('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsMaxGradU.mat','maxGradU','maxGradUTri');

%%


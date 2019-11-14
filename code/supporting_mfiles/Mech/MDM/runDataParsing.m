% clc; clear; close all;
% Folder Directory
% F:\Neutrophil3D\Figures\fig04
clear
%path = 'F:\Neutrophil3D\Data\';
% path = '\\THEYETI\Neutrophil3D\Neutrophil3D\Data\';
path = 'E:\Neutrophil3D\Data';

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
folderDir{16} = 'chemokinesis\Naive\neutr_kineses_May8cell3001';

%%
for i = 1:length(folderDir)
    tic
    cd(fullfile(path,folderDir{i},'Cell Data'));
    disp([num2str(i),'  ', folderDir{i}]);
    
    um2vxl = [0.20716, 0.20716, 0.3];
    dm = 4*um2vxl;
    
    fileDir = dir('Cell*.mat');
    fileDir = {fileDir.name};
    f = cell(length(fileDir),1);
    v = cell(length(fileDir),1);
    
    for j = 1:length(fileDir)
        load(fileDir{j},'FV','vertices');
        f{j} = FV.faces;
        v{j} = bsxfun(@times, vertices, um2vxl);
    end
    
    cd(fullfile(path,folderDir{i}));
    save(fullfile(path,folderDir{i},'resultsCell.mat'), 'f','v')
%     load('resultsDVC.mat');
%     
%     um2vxl = [0.20716, 0.20716, 0.3];
%     dm = 4*um2vxl;
%     save('resultsDVC.mat','u','beadChannel','um2vxl','dm');
%     
%     
%         if exist('MresultsDVC.mat') == 2
%             load('MresultsDVC.mat');
%         elseif exist('resultsDVCgreen.mat') == 2
%             load('resultsDVCgreen.mat');
%         else
%             load('resultsDVC.mat');
%         end
%     
%         dm = 4*um2vxl;
%         savepath = fullfile('E:\Neutrophil3D\Data',folderDir{i},'resultsDVC.mat');
%         save(savepath,'beadChannel','u','um2vxl','dm');
    
end

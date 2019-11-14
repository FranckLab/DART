
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

for i = 4%:length(folderDir)
    cd(fullfile(path,folderDir{i}));
    disp([num2str(i),'  ', folderDir{i}]);
    
    
    load('resultsDVC.mat','u','um2vxl','dm');
    u = changeUnits(u,um2vxl);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('resultsCell.mat');
    nTime   = min([length(u), length(f) - 1]);       % number of DVC increments points
    
    v0      = v{1}; 
    v1      = v; % from stain
    
    uTri    = cell(nTime,1);
    uTri1   = cell(nTime,1);
    v       = cell(nTime + 1, 1);
    
    v{1} = v0; % from DVC data
    for j = 1:nTime
        uTri{j} = triInterp(v{j},u{j},dm);
        v{j + 1} = v{j} + uTri{j};  % move surface with displacements
        uTri1{j} = triInterp(v1{j},u{j},dm);
    end
    
%%

lims = [min(cell2mat([v;v1])) ; max(cell2mat([v;v1]))];
lims = lims(:)';

    close all;
    figure;

    
    for j = 1:nTime
        h = trisurf(f{j},v1{j}(:,1),v1{j}(:,2),v1{j}(:,3),vecMag(uTri1{j},2)); % deformed cell
        set(h,'linestyle','none'); c = colorbar('location','eastoutside');
        ylabel(c,'|u| (\mum)','fontsize',20)
        light; view(3); axis image; axis(lims)
        xlabel('x (\mum)'); ylabel('y (\mum)'); zlabel('z (\mum)');
        grid on;
        set(gca,'color',[0.5,0.5,0.5]);
        title(['Time: ',num2str(j*2),' min [Stain]']);
        drawnow update
        print(gcf,'-dpng ','-opengl',['Stain', num2str(j,'%2.2i')])
    end
    %%
    figure;
    for j = 1:nTime
        h = trisurf(f{1},v{j}(:,1),v{j}(:,2),v{j}(:,3),vecMag(uTri{j},2));
        set(h,'linestyle','none'); c = colorbar('location','eastoutside');
        ylabel(c,'|u| (\mum)','fontsize',20)
        light; view(3); axis image; axis(lims)
        xlabel('x (\mum)'); ylabel('y (\mum)'); zlabel('z (\mum)');
        grid on;
        set(gca,'color',[0.5,0.5,0.5]);
        title(['Time: ',num2str(j*2),' min [DVC]']);
        drawnow
        print(gcf,'-dpng ','-opengl',['DVC', num2str(j,'%2.2i')])
    end
    %%
    
    
lims = [min(cell2mat([v;v1])) ; max(cell2mat([v;v1]))];
lims = lims(:)';
    close all;
    figure;
    
    
    for j = 1:nTime
        clf
        title(['Time: ',num2str(j*2),' min [RED - Stain, BLUE - DVC]']);
        set(gca,'color',[0.5,0.5,0.5]);
        hold on
        h = trisurf(f{j},v1{j}(:,1),v1{j}(:,2),v1{j}(:,3),'facecolor','r','facealpha', 0.5);
        set(h,'linestyle','none');
        
        h = trisurf(f{1},v{j}(:,1),v{j}(:,2),v{j}(:,3),'facecolor','b','facealpha', 0.5);
        set(h,'linestyle','none');
        
        light; view(3); 
        axis image;
        axis(lims)
        xlabel('x (\mum)'); ylabel('y (\mum)'); zlabel('z (\mum)');
        drawnow
        hold off;
        grid on;
        pause(0.01)
        
        
         print(gcf,'-dpng ','-opengl',['Overlay', num2str(j,'%2.2i')])
%         
        
    end
    
    
end
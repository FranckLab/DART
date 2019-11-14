
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

% %% compute ellipse
% for i = 1:length(folderDir)
%     cd(fullfile(path,folderDir{i}));
%     disp([num2str(i),'  ', folderDir{i}]);
%
%     load('resultsCell.mat');
%
%     nTime       =  length(f);
%     fE          =  cell(nTime,1);
%     vE          =  cell(nTime,1);
%     ellipseInfo =  cell(nTime,1);
%
%     for j = 1:nTime
%         [fE{j}, vE{j}, ellipseInfo{j}] = triMinEllipse(v{j},0.30);
%     end
%     save('resultsEllipse.mat','vE','fE','ellipseInfo')
% end

mechanicsE = cell(length(folderDir),1);
mechanics  = cell(length(folderDir),1);
%%

for i = 1:length(folderDir)
    cd(fullfile(path,folderDir{i}));
    disp([num2str(i),'  ', folderDir{i}]);
    
    
    load('resultsDVC.mat','u','um2vxl','dm');
    u0 = u;
    
    stretch = (1:0.05:1.5)';
    nStretch = length(stretch);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate for dilated ellipse
    load('resultsCell.mat');
    nTime   = min([length(u), length(f) - 1]);       % number of DVC increments points
    
    v0 = v;
    
    mechanics = cell(length(stretch),1);
    
    for ii = 1:nStretch
        
        
        uTri    = cell(nTime,1);
        v       = cell(nTime,1);
        
        for j = 1:nTime
            for k = 1:3, u{j}{k} = u0{j}{k}*um2vxl(k); end
             u{j} = removeNoise_Neutrophil3D(u{j}(1:3),2);
             
             
            v{j} = triScale(v0{j},stretch(ii));
            fp = triCenters(f{j},v{j});
            uTri{j}           = triInterp(fp,u{j},dm);
        end
        
        mechanics{ii,1} = calculateMeanMetrics(f,v,uTri);
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate for dilated ellipse
    load('resultsEllipse.mat');
    v0 = vE;
    f = fE;
    
    mechanicsE = cell(length(stretch),1);
    
    for ii = 1:nStretch
        
        uTri    = cell(nTime,1);
        
        for j = 1:nTime
            for k = 1:3, u{j}{k} = u0{j}{k}*um2vxl(k); end
            u{j} = removeNoise_Neutrophil3D(u{j}(1:3),2);
            
            v{j} = triScale(v0{j},stretch(ii));
            fp = triCenters(f{j},v{j});
            uTri{j}           = triInterp(fp,u{j}(:,1:3),dm);
        end
        
        mechanicsE{ii,1} = calculateMeanMetrics(f,v,uTri);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    save('resultsMechanicsBoundary.mat','mechanicsE','mechanics','stretch');
    
end



%% bin data
clearvars -except folderDir path
% nStretch = 11;

% E = cell(nStretch,1);
% EE = cell(nStretch,1);
e = cell(1,16);
eE = cell(1,16);

for i = 1:length(folderDir)
    cd(fullfile(path,folderDir{i}));
    disp([num2str(i),'  ', folderDir{i}]);
    
    
    
    load('resultsMechanicsBoundary.mat');
    load('resultsMechanics.mat');
    nStretch = length(stretch);
    
    E0 = mechanics{stretch == 1}.E;
    V0 = mechanics{stretch == 1}.triV;
    nTime = size(E0,3);
    type = 'fro';
    
    for ii = 1:nStretch
        E = mechanics{ii}.E;
        V = mechanics{ii}.triV;
        
        EE = mechanicsE{ii}.E;
        VE = mechanicsE{ii}.triV;
        
        
        
        for j = 1:nTime
            I = eye(3);
%             F0_  = E0(:,:,j);
%             E0_ = 1/2*( F0_' * F0_ - I);
%             
%             dF_ = E(:,:,j) - F0_ + I;
%             dE_ = 1/2*( dF_' * dF_ - I);
%             
%             dFE_ = EE(:,:,j) - F0_ + I;
%             dEE_ = 1/2*( dFE_' * dFE_ - I);
%             
%             e{ii}(end + 1,1) = norm(dE_ ,type) /  norm(E0_,type);
%              eE{ii}(end + 1,1) = norm(dEE_ ,type) /  norm(E0_,type);
            
              E0_ = E0(:,:,j);
            E_ = E(:,:,j);
            EE_ = EE(:,:,j);
            
                       
%             E0_ = V0(j)* (E0(:,:,j) - I);
%             E_ = V(j) * (E(:,:,j) - I);
%             EE_ = VE(j) * (EE(:,:,j) - I);
            
            
            
            e{ii}(end + 1,1) = norm(E_ - E0_ ,type) /  norm(E0_,type);
             eE{ii}(end + 1,1) = norm(EE_ - E0_ ,type) /  norm(E0_,type);
        end
    end
    
    
%     
%     for ii = 1:nStretch
%         E{ii} = cat(3, E{ii}, mechanics{ii}.E);
%         EE{ii} = cat(3, EE{ii}, mechanicsE{ii}.E);
%     end
%     
    
end

e = cell2mat(e);
eE = cell2mat(eE);

% prctile(e,25)
prctile(e,50)

prctile(eE,50)

%%
idx = stretch >= 1;
LB = 25;
UB = 75;
figure; 
hold on;
% plot(stretch(idx),prctile(e(:,idx),LB),'-r',stretch(idx),prctile(eE(:,idx),LB),'-b');
% plot(stretch(idx),prctile(e(:,idx),UB),'-r',stretch(idx),prctile(eE(:,idx),UB),'-b');
jbfill(stretch(idx)',prctile(eE(:,idx),UB),prctile(eE(:,idx),LB),[0 0 255]/255,'none',0,0.5)
jbfill(stretch(idx)',prctile(e(:,idx),UB),prctile(e(:,idx),LB),[85 85 85]/255,'none',0,1)
plot(stretch(idx),prctile(e(:,idx),50),'-k','linewidth',1.5);
plot(stretch(idx),prctile(eE(:,idx),50),':k','linewidth',1.5);
hold off


set(gca, 'xtick', [1:0.125:1.5]);

set(gca, 'ytick',[0:0.25:1]);
set(gca, 'Ticklength',[0 0])


h = gcf;
set(h, 'PaperUnits','inches')
set(h, 'PaperSize',[3 2])
set(h, 'Paperposition', [3 3 3 2]);
set(h, 'Units','inches','Position', [3, 3, 3 2]);

%
% savepath = 'C:\Users\Eyal\Google Drive\Research\Publications\14 Neutrophil3D\manuscript\figures\SI\boundaryError';
% savename = fullfile(savepath,'boundaryError_plot.eps');
% export_fig(h, savename,  '-eps', '-painters');
% print(h,   '-deps', '-painters', savename);

% prctile(e,75)

%% calculate error

% E0 = E{stretch == 1};









% %%
%     nStretch = length(stretch);
%     nTime = length(mechanics{1}.filename);
%     eType = 'fro';
%
%     eData0 = mechanics{1}.E;
%     eFro = zeros(nTime,nStretch);
%     eFroE = zeros(nTime,nStretch);
%
%     for l = 1:nStretch
%
%         eData = mechanics{l}.E;
%         eDataE = mechanicsE{l}.E;
%
%         for k = 1:nTime
%
%             eFro(k,l) = norm(eData(:,:,k) - eData0(:,:,k),eType) / norm(eData0(:,:,k),eType);
%             eFroE(k,l) = norm(eDataE(:,:,k) - eData0(:,:,k),eType) / norm(eData0(:,:,k),eType);
% %             eFro(k,l) = norm(eData(k,:) - eData0(k,:),eType);
% %             eFroE(k,l) = norm(eDataE(k,:) - eData0(k,:),eType);
% %             eFro(k,l) = norm(eData(k,:) - eData0(k,:),eType) / norm(eData0(k,:),eType);
% %             eFroE(k,l) = norm(eDataE(k,:) - eData0(k,:),eType) / norm(eData0(k,:),eType);
%         end
%     end
% 
%
%     subplot(4,4,i)
% hold on;
% errorbar(stretch,prctile(eFro,50),prctile(eFro,25),prctile(eFro,75),'r.-')
% errorbar(stretch, prctile(eFroE,50),prctile(eFroE,25),prctile(eFroE,75),'b*-')
% hold off;
%
% fontSize = 8;
%
%
% set(gca,'fontsize',fontSize,'xtick',[0.5:0.5:3.5])
% hl = legend('cell Dilation','ellipse Dilation');
% set(hl,'FontSize',fontSize)
% xlabel('stretch','fontsize',fontSize);
% ylabel('median( ||<Ei> - <E''i>||_\infty/||<Ei>||_\infty )','fontsize',fontSize)
% title(folderDir{i},'fontsize',fontSize);
% box on;
%
% end
%
%
%
%
%
% % cd(fullfile(path,folderDir{1},'raw'));
% % dirI = dir('dec*.mat');
% % dirI = {dirI.name};
% %
% % for i = 1:length(dirI)
% %    load(dirI{i});
% %
% %    I = double(vol{1});
% %    I = I/max(I(:));
% %    um2vxl = [0.21,0.21,0.3];
% %    figure; image2(max(vol{1}(:,:,round(lim(1):lim(2))),[],3))
% %    figure; image2(max(vol{2}(:,:,round(lim(1):lim(2))),[],3))
% % end

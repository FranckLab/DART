
% Folder Directory

clc; clear; close all;
% path = 'F:\Neutrophil3D\Data\';
% path = '\\THEYETI\Neutrophil3D\Neutrophil3D\Data\';
% path = 'E:\Neutrophil3D\Data';
% path = 'C:\Users\Eyal\Dropbox\Neutrophil3D\data';
path = 'E:\Jon\Google Drive\14 Neutrophil3D\data';
%folderDir{01} = 'mFiles\validateEshelby'

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

%% compute sphere
for i = 1:length(folderDir)
    cd(fullfile(path,folderDir{i}));
    disp([num2str(i),'  ', folderDir{i}]);

    load('resultsCell.mat');
    load('resultsEllipse.mat');
    nTime       =  length(f);
    fS          =  cell(nTime,1);
    vS          =  cell(nTime,1);
    sphereInfo =  cell(nTime,1);

    for j = 1:nTime
        %[fS{j}, vS{j}, sphereInfo{j}] = triMinEllipse(v{j},0.30);
        sphereInfo{j} = ellipseInfo{j};
        sphereInfo{j}.radius = max(ellipseInfo{j}.radius);
        [fS{j}, vS{j}] = sphere_tri('ico',4,sphereInfo{j}.radius);
    end
    save('resultsSphere.mat','vS','fS','sphereInfo')
end

%%

for i = 1:length(folderDir)
    cd(fullfile(path,folderDir{i}));
    disp([num2str(i),'  ', folderDir{i}]);
    
    
    load('resultsDVC.mat','u','um2vxl','dm');
    %load('resultsDVC_Esh.mat','u','um2vxl','dm');
    u0 = u;
    
    %stretch = (1:0.1:3)';
    %stretch = (0.5:0.1:2.5)';
    stretch = 1;
    nStretch = length(stretch);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate for dilated ellipse
    load('resultsCell.mat');
    %load('resultsCell_Esh.mat'); a{1} = v; v=a; a{1} = f; f=a;
    
    nTime   = min([length(u), length(f) - 1]);       % number of DVC increments points
    %nTime = 1; %Added for Esh solution
    
    v0 = v;
    
    mechanics = cell(length(stretch),1);
    
    for l = 1:nStretch
        
        uTri    = cell(nTime,1);
        v       = cell(nTime,1);
        
        for j = 1:nTime
            for k = 1:3, u{j}{k} = u0{j}{k}*um2vxl(k); end
            
            v{j} = triScale(v0{j},stretch(l));
            fp = triCenters(f{j},v{j});
            uTri{j}           = triInterp(fp,u{j}(:,1:3),dm);
        end
        
        mechanics{l} = calculateMeanMetrics(f,v,uTri);
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate for dilated ellipse
    load('resultsEllipse.mat');
    v0 = vE;
    f = fE;
    
    mechanicsE = cell(length(stretch),1);
    
    for l = 1:nStretch
        
        uTri    = cell(nTime,1);
        
        for j = 1:nTime
            for k = 1:3, u{j}{k} = u0{j}{k}*um2vxl(k); end
            
            v{j} = triScale(v0{j},stretch(l));
            fp = triCenters(f{j},v{j});
            uTri{j}           = triInterp(fp,u{j}(:,1:3),dm);
        end
        
        mechanicsE{l} = calculateMeanMetrics(f,v,uTri);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate for dilated ellipse
    load('resultsSphere.mat');
    v0 = vS;
    f = fS;
    
    mechanicsS = cell(length(stretch),1);
    
    for l = 1:nStretch
        
        uTri    = cell(nTime,1);
        
        for j = 1:nTime
            for k = 1:3, u{j}{k} = u0{j}{k}*um2vxl(k); end
            
            v{j} = triScale(v0{j},stretch(l));
            fp = triCenters(f{j},v{j});
            uTri{j}           = triInterp(fp,u{j}(:,1:3),dm);
        end
        
        mechanicsS{l} = calculateMeanMetrics(f,v,uTri);
    end
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save('resultsMechanicsNES.mat','mechanicsS','mechanicsE','mechanics','stretch');
    
end



%% calculate error

estorFro = [];%zeros(length(folderDir),1);
estorFroE = [];%zeros(length(folderDir),1);
estorFroS = [];%zeros(length(folderDir),1);

figure;
for i = 1:length(folderDir)
    cd(fullfile(path,folderDir{i}));
    disp([num2str(i),'  ', folderDir{i}]);
    
    load('resultsMechanicsNES.mat');
    nStretch = length(stretch);
    nTime = length(mechanics{1}.triV);%length(mechanics{1}.filename);
    eType = 'fro';
    
    %eData0 = mechanics{1}.F;
    eFro = zeros(nTime,nStretch);
    eFroE = zeros(nTime,nStretch);
    eFroS = zeros(nTime,nStretch);
    
    for l = 1:nStretch
        
        eData = mechanics{l}.F;
        eDataE = mechanicsE{l}.F;
        eDataS = mechanicsS{l}.F;
        
        for k = 1:nTime
            %eFro(k,l) = norm(eData(:,:,k),eType);
            %eFroE(k,l) = norm(eDataE(:,:,k),eType);
            %eFroS(k,l) = norm(eDataS(:,:,k),eType);
            eFroE(k,l) = norm(eDataE(:,:,k) - eData(:,:,k),eType) / norm(eData(:,:,k),eType);
            eFroS(k,1) = norm(eDataS(:,:,k) - eData(:,:,k),eType) / norm(eData(:,:,k),eType);
            
            %eFro(k,l) = norm(eData(:,:,k) - eData0(:,:,k),eType) / norm(eData0(:,:,k),eType);
            %eFroE(k,l) = norm(eDataE(:,:,k) - eData0(:,:,k),eType) / norm(eData0(:,:,k),eType);
            % eFro(k,l) = norm(eData(:,:,k) - eData0(:,:,k),eType);
            % eFroE(k,l) = norm(eDataE(:,:,k) - eData0(:,:,k),eType);
            %eFro(k,l) = norm(eData(:,:,k)- eye(3,3),eType); %using grad u 
            %eFroE(k,l) = norm(eDataE(:,:,k)- eye(3,3),eType);
            %             eFro(k,l) = norm(eData(k,:) - eData0(k,:),eType) / norm(eData0(k,:),eType);
            %             eFroE(k,l) = norm(eDataE(k,:) - eData0(k,:),eType) / norm(eData0(k,:),eType);
        end
    end
    
    
    subplot(4,4,i)
    hold on;
    errorbar(stretch,prctile(eFro,50),prctile(eFro,25),prctile(eFro,75),'r.-')
    errorbar(stretch,prctile(eFroE,50),prctile(eFroE,25),prctile(eFroE,75),'b*-')
    hold off;
    estorFro = [estorFro; eFro];
    estorFroE = [estorFroE; eFroE];%prctile(eFro,50);
    estorFroS = [estorFroS; eFroS];%prctile(eFroE,50);
    
    fontSize = 8;
    
    
    set(gca,'fontsize',fontSize,'xtick',[0.5:0.5:3.5])
    hl = legend('cell Dilation','ellipse Dilation');
    set(hl,'FontSize',fontSize)
    xlabel('stretch','fontsize',fontSize);
    ylabel('median( ||<Fi> - <F''i>||_\infty/||<Fi>||_\infty )','fontsize',fontSize)
    %ylabel('median( ||<gradU>||)','fontsize',fontSize);
    title(folderDir{i},'fontsize',fontSize);
    box on;
    
end

figure;
    hold on;
    errorbar(1,prctile(estorFroE,50),prctile(estorFroE,25),prctile(estorFroE,75),'r.-')
    errorbar(2,prctile(estorFroS,50),prctile(estorFroS,25),prctile(estorFroS,75),'b*-')
    xlim([0,3]);
    hold off;

%%
figure;
plot(stretch,median(estorFroE),'--','Color','b')
hold on;
plot(stretch,prctile(estorFroE,25),'-','Color','b')
plot(stretch,prctile(estorFroE,75),'-','Color','b')
set(gca,'fontsize',fontSize,'xtick',[0.5:0.5:3.5])
xlabel('stretch','fontsize',fontSize);
%ylabel('median( ||<Ei>||)','fontsize',fontSize);
ylabel('median( ||<Fi> - <F''i>||_\infty/||<Fi>||_\infty )','fontsize',fontSize)
title('All Cells','fontsize',fontSize);

figure;
plot(stretch,median(estorFro),'--','Color','r')
hold on;
plot(stretch,prctile(estorFro,25),'-','Color','r')
plot(stretch,prctile(estorFro,75),'-','Color','r')
set(gca,'fontsize',fontSize,'xtick',[0.5:0.5:3.5])
xlabel('stretch','fontsize',fontSize);
%ylabel('median( ||<Ei>||)','fontsize',fontSize);
ylabel('median( ||<Fi> - <F''i>||_\infty/||<Fi>||_\infty )','fontsize',fontSize)
title('All Cells','fontsize',fontSize);

% cd(fullfile(path,folderDir{1},'raw'));
% dirI = dir('dec*.mat');
% dirI = {dirI.name};
%
% for i = 1:length(dirI)
%    load(dirI{i});
%
%    I = double(vol{1});
%    I = I/max(I(:));
%    um2vxl = [0.21,0.21,0.3];
%    figure; image2(max(vol{1}(:,:,round(lim(1):lim(2))),[],3))
%    figure; image2(max(vol{2}(:,:,round(lim(1):lim(2))),[],3))
% end

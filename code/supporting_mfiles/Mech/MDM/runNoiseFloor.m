
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
 u1 = [];
    u2 = [];
    u3 = [];
    
for i = 1:length(folderDir)
    cd(fullfile(path,folderDir{i}));
    disp([num2str(i),'  ', folderDir{i}]);
    
    
    load('resultsDVC.mat','u','um2vxl','dm');
    u = changeUnits(u,um2vxl);
    nTime = length(u);
    
   
    i
    for j = 1:nTime
        for k = 1:3
            s = 4;
            u{j}{k} = u{j}{k}(4: (end - 4), 4: (end - 4), 4: (end - 4));
        end
            
        u1 = [u1; u{j}{1}(:)];
        u2 = [u2; u{j}{2}(:)];
        u3 = [u3; u{j}{3}(:)];
        
        
        
        
%         uSTD = std([u1 u2 u3]);
%         uCount = numel(u1);
        
        
        
        
        
%         [hu,pvalue] = kstest((u1 - mean(u1))/std(u1);
        
        
        
        
    end
        
    
    
    
%     
%     I{i} = I0{i} - mean(noise);
%     noiseLevel{i} = 2*std(noise);
%     
%     I{i}(abs(I{i})<abs(noiseLevel{i})) = 0;
    
    
    
    
end

%% create plots
% Plot results for large sample size in same figure window as the
% small-size distribution of the same type
ksText ={'normal', 'not normal');

for i = 1:length(data);
figure;
Data = data{i};
STD = std(Data);
Mean = mean(Data);
Median = median(Data);
histogram(Data,1000);
% [h, p] = kstest(Data);

axis([(-5*STD + Mean) (5*STD + Mean), get(gca,'YLim')])



title(sprintf('Normally distributed data, %d samples',length(Data)));
text(0,.8,sprintf('Mean = %.1f\nStd. Dev. = %.1f\nMedian = %.1f\n', ...
    Mean, STD, Median));
%         text(0,0.3, sprintf('KS: p = %.3f  %s\n%s',...
%            h,p);

end



%%

figure; 
subplot(3,1,1);

STD = std(u1);
Mean = mean(u1);
hist(u1,10000);
axis([(-5*STD + Mean) (5*STD + Mean), get(gca,'YLim')])
title('histogram of u1');

subplot(3,1,2);
STD = std(u2);
Mean = mean(u2);
hist(u2,10000);
axis([(-5*STD + Mean) (5*STD + Mean), get(gca,'YLim')])
title('histogram of u2');

subplot(3,1,3);
STD = std(u3);
Mean = mean(u3);
hist(u3,10000);
axis([(-5*STD + Mean) (5*STD + Mean), get(gca,'YLim')])
title('histogram of u3');








clc
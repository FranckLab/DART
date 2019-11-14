%% <Ei/|E|> )
clear
load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsMechanics.mat');

histVar = cell(1,length(mechanics));
for i = 1:length(mechanics)
    histVar{i} = (mechanics{i}.eigValENorm);
end

histData = cell(1,4);
histData{1} =     histVar(12:15)'; % Kinesis Naive
histData{2} =     histVar(4:6)';   % taxis Naive
histData{3} =     histVar(7:11)';  % Kinesis LPS
histData{4} =     histVar(1:3)';   % taxis LPS


histData = cellfun(@cell2mat, histData,'uniformoutput',0);

boxData = nan(max(cellfun(@(x) size(x,1), histData)),12);

k = 1;
for j = 1:3
    for i = 1:4
        boxData(1:length(histData{i}(:,j)),k) = histData{i}(:,j);
        k = k + 1;
    end
end

close all;
figure;
boxplot(boxData,'symbol','');
set(gca,'ylim',[-1 1],'ytick',linspace(-1,1,5));
formatBoxPlot([1.4, 1]);
clc

saveDir = 'C:\Users\Eyal\Google Drive\Research\Publications\14 Neutrophil3D\manuscript\figures\boxAndWhiskers\';
export_fig(gcf,[saveDir, 'EiNorm'],'-pdf','-q100','-nocrop');

%% <Ei>
clear
load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsMechanicsNoNoise.mat');

histVar = cell(1,length(mechanics));
for i = 1:length(mechanics)
    histVar{i} = (mechanics{i}.eigValE);
end

histData = cell(1,4);
histData{1} =     histVar(12:15)'; % Kinesis Naive
histData{2} =     histVar(4:6)';   % taxis Naive
histData{3} =     histVar(7:11)';  % Kinesis LPS
histData{4} =     histVar(1:3)';   % taxis LPS

histData = cellfun(@cell2mat, histData,'uniformoutput',0);

boxData = nan(max(cellfun(@(x) size(x,1), histData)),12);

k = 1;
for j = 1:3
    for i = 1:4
        boxData(1:length(histData{i}(:,j)),k) = histData{i}(:,j);
        k = k + 1;
    end
end


close all;
figure;
boxplot(boxData,'symbol','');
set(gca,'ylim',[-0.1 0.2],'ytick',linspace(-0.1,0.2,4));
formatBoxPlot([1.4, 1]);
clc
%
saveDir = 'C:\Users\Eyal\Google Drive\Research\Publications\14 Neutrophil3D\manuscript\figures\boxAndWhiskers\';
export_fig(gcf,[saveDir, 'Ei'],'-pdf','-q100','-nocrop');

%% <J>
clear
% load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsMechanics.mat');
load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsMechanicsNoNoise.mat');

histVar = cell(1,length(mechanics));
for i = 1:length(mechanics)
    histVar{i} = (mechanics{i}.J);
end

histData = cell(1,4);
histData{1} =     histVar(12:15)'; % Kinesis Naive
histData{2} =     histVar(4:6)';   % taxis Naive
histData{3} =     histVar(7:11)';  % Kinesis LPS
histData{4} =     histVar(1:3)';   % taxis LPS

histData = cellfun(@cell2mat, histData,'uniformoutput',0);

boxData = nan(max(cellfun(@(x) size(x,1), histData)),4);

k = 1;
for i = 1:4
    boxData(1:length(histData{i}),k) = histData{i};
    k = k + 1;
end

% boxColor = repmat('k',4,1);

h = figure;
boxplot(boxData,'symbol','', 'widths', 0.5);
set(gca,'ylim',[0.8 1.20],'ytick',[0.80 .9 1 1.1 1.2]);
formatBoxPlot(h, [0.8, 1]);


clc
%
saveDir = 'C:\Users\Eyal\Google Drive\Research\Publications\14 Neutrophil3D\manuscript\figures\boxAndWhiskers\';
export_fig(gcf,[saveDir, 'J'],'-pdf','-q100','-nocrop');
% print(gcf,'-dpdf','-painters',[saveDir, 'J']);

%% <gradU>

clear
load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsMaxGradU.mat');
% load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsMechanicsNoNoise.mat');

histVar = cell(1,length(maxGradU));
for i = 1:length(maxGradU)
    histVar{i} = maxGradUTri{i};
end

histData = cell(1,4);
histData{1} =     histVar(12:15)'; % Kinesis Naive
histData{2} =     histVar(4:6)';   % taxis Naive
histData{3} =     histVar(7:11)';  % Kinesis LPS
histData{4} =     histVar(1:3)';   % taxis LPS

histData = cellfun(@cell2mat, histData,'uniformoutput',0);

boxData = nan(max(cellfun(@(x) size(x,1), histData)),4);

k = 1;
for i = 1:4
    boxData(1:length(histData{i}),k) = histData{i};
    k = k + 1;
end

figure;
boxplot(boxData,'symbol','');
set(gca,'ylim',[0 0.8],'ytick',linspace(0,0.8,5));
formatBoxPlot([0.8, 1]);
clc

saveDir = 'C:\Users\Eyal\Google Drive\Research\Publications\14 Neutrophil3D\manuscript\figures\boxAndWhiskers\';
export_fig(gcf,[saveDir, 'maxGradU'],'-pdf','-q100','-nocrop');

%% <theta>
clear
% load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsMechanics.mat');
load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsMechanicsNoNoise.mat');

histVar = cell(1,length(mechanics));
for i = 1:length(mechanics)
    histVar{i} = (mechanics{i}.axisAngleR(:,4))*180/pi;
end

histData = cell(1,4);
histData{1} =     histVar(12:15)'; % Kinesis Naive
histData{2} =     histVar(4:6)';   % taxis Naive
histData{3} =     histVar(7:11)';  % Kinesis LPS
histData{4} =     histVar(1:3)';   % taxis LPS

histData = cellfun(@cell2mat, histData,'uniformoutput',0);

boxData = nan(max(cellfun(@(x) size(x,1), histData)),4);

k = 1;
for i = 1:4
    boxData(1:length(histData{i}),k) = histData{i};
    k = k + 1;
end

close all;
boxColor = repmat('k',4,1);

figure;
boxplot(boxData,'symbol','');
set(gca,'ylim',[0 3.2],'ytick',linspace(0,3.2,5));
formatBoxPlot([0.8, 1]);
clc

saveDir = 'C:\Users\Eyal\Google Drive\Research\Publications\14 Neutrophil3D\manuscript\figures\boxAndWhiskers\';
export_fig(gcf,[saveDir, 'theta'],'-pdf','-q100','-nocrop');
% print(gcf,'-dpdf','-painters',[saveDir, 'theta']);

%% <int(theta)>
clear
% load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsMechanics.mat');
load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsMechanicsNoNoise.mat');

dt = 2;
histVar = cell(1,length(mechanics));
for i = 1:length(mechanics)
    histVar{i} = (mechanics{i}.axisAngleR(:,4))*180/pi;
end

histData = cell(1,4);
histData{1} =     histVar(12:15)'; % Kinesis Naive
histData{2} =     histVar(4:6)';   % taxis Naive
histData{3} =     histVar(7:11)';  % Kinesis LPS
histData{4} =     histVar(1:3)';   % taxis LPS


intData = cell(1,4);
for i = 1:length(histData)
    for j = 1:length(histData{i})
        
        t = (0:dt:dt*length(histData{i}{j}(:,1)));
        theta = [0; histData{i}{j}(:,1)];
        intData{i}{j}(:,1) = t;
        intData{i}{j}(:,2) = cumtrapz(t,abs(theta)); % calculate cumalitive sum
        
    end
end

linespec = {'-',':','-',':'};
linecolor = {'r','r','b','b'};



for i = 1:4
    intSize = sort(cellfun(@(x) size(x,1), intData{i}));
    intSize = intSize(end - 1);
    
    
    
    int = nan(intSize, size(intData{i},2));
    
    for ii = 1:length(intData{i})
        if size(intData{i}{ii},1) <= intSize
             idx = 1:size(intData{i}{ii},1);
        else
            idx = 1:intSize;
        end
        int(idx,ii) = intData{i}{ii}(idx,2);
    end
    
    time = (0:dt:dt*(size(int,1))-1)';
    intMean = nanmean(int,2);
    intStd = nanstd(int,1,2);
    intLines{i} = [time, intMean - intStd, intMean, intMean + intStd];
    
end

linecolor = {'r','b','m','c'};

figure;
for i = 1:4
%     for ii = 2:size(intLines{i},2)
        hold on;
        plot(intLines{i}(:,1),intLines{i}(:,3),'linewidth',2,'color',linecolor{i});
        jbfill(intLines{i}(:,1)', intLines{i}(:,4)', intLines{i}(:,2)',linecolor{i})
%         plot(intLines{i}(:,1),intLines{i}(:,ii),'linewidth',2,'color',linecolor{i});
        hold off;
%     end

end


saveDir = 'C:\Users\Eyal\Google Drive\Research\Publications\14 Neutrophil3D\manuscript\figures\boxAndWhiskers\';
print(gcf,'-dpdf','-painters',[saveDir, 'intTheta']);
% xlabel('time (min)');
% ylabel('\int_0^t|\theta(\tau)|d\tau (\circ)');
% 
%Add legend
% h=legend('T - N', 'K - N', 'T - L', 'K - L');

%% <|v|>
% clear
% % load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsCentroid.mat');
% load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsMechanicsNoNoise.mat');
% 
% histVar = cell(1,length(centroid));
% for i = 1:length(centroid)
%     C = centroid{i}.CA;
%     dCdt = diff(C)/2;
%     histVar{i} = vecMag(dCdt,2);
%     
% end
% 
% histData = cell(1,4);
% histData{1} =     histVar(12:15)'; % Kinesis Naive
% histData{2} =     histVar(4:6)';   % taxis Naive
% histData{3} =     histVar(7:11)';  % Kinesis LPS
% histData{4} =     histVar(1:3)';   % taxis LPS
% 
% 
% histData = cellfun(@cell2mat, histData,'uniformoutput',0);
% 
% sizeHistData = reshape(cell2mat(cellfun(@size, histData,'uniformoutput',0)),2,[]);
% boxData = nan(max(sizeHistData(1,:)), sum(sizeHistData(2,:)));
% 
% k = 1;
% for j = sizeHistData(2)
%     for i = 1:size(histData,2)
%         boxData(1:length(histData{i}(:,j)),k) = histData{i}(:,j);
%         k = k + 1;
%     end
% end
% 
% boxColor = repmat('k',4,1);
% 
% close all;
% figure;
% boxplot(boxData,'symbol','','color',boxColor);
% set(gca,'ylim',2.1*[0 1],'ytick',[0,1 2]);
% formatBoxPlot();
% clc
% 
% saveDir = 'C:\Users\Eyal\Google Drive\Research\Publications\14 Neutrophil3D\manuscript\figures\boxAndWhiskers\';
% export_fig(gcf,[saveDir, 'speed'],'-pdf','-q100','-nocrop');
% print(gcf,'-dpdf','-painters',[saveDir, 'speed']);

%% phi
% clear
% % load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsCentroid.mat');
% load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsMechanicsNoNoise.mat');
% 
% histVar = cell(1,length(centroid));
% 
% for i = 1:length(centroid)
%     C = centroid{i}.CA;
%     dC = diff(C);
%     cosPhi = (vecCos(dC(1:end-1,:),dC(2:end,:),2));
%     histVar{i} = cosPhi;
%     
% end
% 
% histData = cell(1,4);
% histData{1} =     histVar(12:15)'; % Kinesis Naive
% histData{2} =     histVar(4:6)';   % taxis Naive
% histData{3} =     histVar(7:11)';  % Kinesis LPS
% histData{4} =     histVar(1:3)';   % taxis LPS
% 
% 
% 
% histData = cellfun(@cell2mat, histData,'uniformoutput',0);
% 
% sizeHistData = reshape(cell2mat(cellfun(@size, histData,'uniformoutput',0)),2,[]);
% boxData = nan(max(sizeHistData(1,:)), sum(sizeHistData(2,:)));
% 
% k = 1;
% for j = sizeHistData(2)
%     for i = 1:size(histData,2)
%         boxData(1:length(histData{i}(:,j)),k) = acosd(histData{i}(:,j));
%         k = k + 1;
%     end
% end
% 
% % close all;
% boxColor = repmat('k',4,1);
% 
% figure;
% boxplot(boxData,'symbol','','color',boxColor);
% % set(gca,'ylim',0.45*[0 1],'ytick',[0,0.2,0.4]);
% formatBoxPlot();
% clc
% 
% % saveDir = 'C:\Users\Eyal\Google Drive\Research\Publications\14 Neutrophil3D\manuscript\figures\boxAndWhiskers\';
% % export_fig(gcf,[saveDir, 'CI'],'-pdf','-q100','-nocrop');
% % print(gcf,'-dpdf','-painters',[saveDir, 'phi']);

%% <R> vs <U>
clear
% load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsMechanics.mat');
load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsMechanicsNoNoise.mat');

histVar = cell(1,length(mechanics));
for i = 1:length(mechanics)
    F = mechanics{i}.F;
     
    norme = zeros(size(F,3),1);
    normw = zeros(size(F,3),1);
    normGradU = zeros(size(F,3),1);
    
    type = 'fro';
    I = eye(3,3);
    for ii = 1:size(F,3) 
        w = 1/2*(F(:,:,ii) - F(:,:,ii)');
        e = 1/2*(F(:,:,ii) + F(:,:,ii)' - 2*I);
        
        normw(ii) = norm(w,type);
        norme(ii) = norm(e,type);
        normGradU(ii) = norm(F(:,:,ii) - I,type);
    end
    
%     normSum = normw + norme;
    normSum = normGradU;
    ratiow = normw./normSum;
    ratioe = norme./normSum;
% ratiow = normw;
% ratioe = norme;
    histVar_w{i} = ratiow;
    histVar_e{i} = ratioe;
%     temp{i} = normGradU
    
%     % R AND U
%     R = mechanics{i}.R;
%     U = mechanics{i}.U;
%     
%     normR = zeros(size(R,3),1);
%     normU = zeros(size(U,3),1);
%     I = eye(3,3);
%     type = 'fro';
%     for ii = 1:size(R,3)
%         normR(ii) = norm(R(:,:,ii) - I,type);
%         normU(ii) = norm(U(:,:,ii) - I,type);
%     end
%     
%     normSum = (normR + normU)
%     ratioR = normR./normSum;
%     ratioU = normU./normSum;
%     histVar{i} = ratioR;
    
    
    
end

histData = cell(1,8);
histData{1} =     histVar_w(12:15)'; % Kinesis Naive
histData{2} =     histVar_w(4:6)';   % taxis Naive
histData{3} =     histVar_w(7:11)';  % Kinesis LPS
histData{4} =     histVar_w(1:3)';   % taxis LPS
histData{5} =     histVar_e(12:15)'; % Kinesis Naive
histData{6} =     histVar_e(4:6)';   % taxis Naive
histData{7} =     histVar_e(7:11)';  % Kinesis LPS
histData{8} =     histVar_e(1:3)';   % taxis LPS

histData = cellfun(@cell2mat, histData,'uniformoutput',0);

boxData = nan(max(cellfun(@(x) size(x,1), histData)),8);

k = 1;
for i = 1:size(histData,2)
    boxData(1:length(histData{i}),k) = histData{i};
    k = k + 1;
end


close all;
figure;
hold on;
boxplot(boxData,'symbol','');
% set(gca,'ylim',[0, 1], 'ytick', linspace(0, 1, 5));
% formatBoxPlot([1.5 1]);

% saveDir = 'C:\Users\Eyal\Google Drive\Research\Publications\14 Neutrophil3D\manuscript\figures\boxAndWhiskers\';
% export_fig(gcf,[saveDir, 'sym_vs_skew'],'-pdf','-q100','-nocrop');

%% eigenvalues of <U>

clear
load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsMechanicsNoNoise.mat');

histVar = cell(1,length(mechanics));
for i = 1:length(mechanics)
    histVar{i} = (mechanics{i}.eigValE);
end

histData = cell(1,4);
histData{1} =     histVar(12:15)'; % Kinesis Naive
histData{2} =     histVar(4:6)';   % taxis Naive
histData{3} =     histVar(7:11)';  % Kinesis LPS
histData{4} =     histVar(1:3)';   % taxis LPS

histData = cellfun(@cell2mat, histData,'uniformoutput',0);

boxData = nan(max(cellfun(@(x) size(x,1), histData)),12);

k = 1;
for j = 1:3
    for i = 1:4
        boxData(1:length(histData{i}(:,j)),k) = psqrt(2*histData{i}(:,j)+1);
        k = k + 1;
    end
end


close all;
figure;
boxplot(boxData,'symbol','');
set(gca,'ylim',[-0.1 0.2]+1,'ytick',linspace(-0.1+1,0.2+1,4));
formatBoxPlot([1.4, 1]);
clc
%
saveDir = 'C:\Users\Eyal\Google Drive\Research\Publications\14 Neutrophil3D\manuscript\figures\boxAndWhiskers\';
export_fig(gcf,[saveDir, 'Ei'],'-pdf','-q100','-nocrop');

%% <gradu>
clear
% load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsMechanics.mat');
load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsMechanicsNoNoise.mat');

histVar = cell(1,length(mechanics));
for i = 1:length(mechanics)
    F = mechanics{i}.F;
     
    normGradU = zeros(size(F,3),1);
    
    type = 'fro';
    I = eye(3,3);
    for ii = 1:size(F,3) 
        normGradU(ii) = norm(F(:,:,ii) - I,type);
    end
 
    histVar{i} = normGradU;
end

histData = cell(1,4);
histData{1} =     histVar(12:15)'; % Kinesis Naive
histData{2} =     histVar(4:6)';   % taxis Naive
histData{3} =     histVar(7:11)';  % Kinesis LPS
histData{4} =     histVar(1:3)';   % taxis LPS

histData = cellfun(@cell2mat, histData,'uniformoutput',0);
boxData = nan(max(cellfun(@(x) size(x,1), histData)),4);

k = 1;
for i = 1:size(histData,2)
    boxData(1:length(histData{i}),k) = histData{i};
    k = k + 1;
end


close all;
figure;
hold on;
boxplot(boxData,'symbol','');
set(gca,'ylim',[0, 0.3], 'ytick', linspace(0, 0.3, 5));
% formatBoxPlot([1.5 1]);

% saveDir = 'C:\Users\Eyal\Google Drive\Research\Publications\14 Neutrophil3D\manuscript\figures\boxAndWhiskers\';
% export_fig(gcf,[saveDir, 'sym_vs_skew'],'-pdf','-q100','-nocrop');


%% OLD
% %% <chemotactic index>
% clear
% load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsCentroid.mat');
% % http://web.mit.edu/course/other/beh.410/www/Handouts/motility.pdf
% % http://onlinelibrary.wiley.com/doi/10.1002/aic.690391210/pdf
% % http://climate.indiana.edu/RobesonPubs/hanson_etal_vector_correlation.pdf
% 
% histVar = cell(1,length(centroid));
% factor = zeros(length(centroid),1);
% CI = zeros(length(centroid),1);
% 
% for i = 1:length(centroid)
%     C = centroid{i}.CA;
%     dC = diff(C);
%     totalDistance = (vecMag(dC,2));
%     
%     phi = acos(vecCos(dC(1:end-1,:),dC(2:end,:),2));
%     
%     totalTime = 2*(numel(phi) - 1);
%     t0 = 0:2:totalTime;
%     clear phiRMS
%     dt = 1e-5;
%     dt_ = dt:dt:dt*1000;
%     for ii = 1:length(dt_)
%         t = 0:dt_(ii):totalTime;
%         phi_ = interp1(t0, phi', t,'spline');
%         phiRMS(ii) = sqrt(1/numel(phi_)*sum(phi_.^2));
%     end
%     
%     
%     phiRMS2 = phiRMS.^2;
%     
%     %         figure; plot( dt_,2./phiRMS2);
%     
%     p(i,:) = polyfit(dt_(2:100),2./phiRMS2(2:100),1);
%     persistenceTime = p(i,2);
%     
%     %     dx = abs(sum(dC(:,1)));
%     %     L = norm(sum(dC));
%     %
%     %         factor = (1 - (2/persistenceTime)^(-1)*(1 - exp(-2/persistenceTime)));
%     %         CI(i,1) = dx/L;
%     %
%     totalTime = 2*(size(C,1)-1);
%     t0 = 0:2:totalTime;
%     t = 0:persistenceTime:totalTime;
%     clear C_;
%     C_(:,1) = interp1(t0, C(:,1)', t,'spline');
%     C_(:,2) = interp1(t0, C(:,2)', t,'spline');
%     C_(:,3) = interp1(t0, C(:,3)', t,'spline');
%     
%     
%     dC_ = diff(C_);
%     %      dx = abs(dC_(:,3));
%     %     L = vecMag(dC_,2);
%     CI = abs(vecNorm(dC_,2));
%     
%     %     CI = dx./L;
%     histVar{i} = CI(:,1);
%     
% end
% 
% histData = cell(1,2);
% histData{1} =     histVar(4:6)';   % taxis Naive
% % histData{2} =     histVar(12:15)'; % Kinesis Naive
% 
% histData{2} =     histVar(1:3)';   % taxis LPS
% % histData{4} =     histVar(7:11)';  % Kinesis LPS
% 
% 
% histData = cellfun(@cell2mat, histData,'uniformoutput',0);
% 
% sizeHistData = reshape(cell2mat(cellfun(@size, histData,'uniformoutput',0)),2,[]);
% boxData = nan(max(sizeHistData(1,:)), sum(sizeHistData(2,:)));
% 
% k = 1;
% for j = sizeHistData(2)
%     for i = 1:size(histData,2)
%         boxData(1:length(histData{i}(:,j)),k) = histData{i}(:,j);
%         k = k + 1;
%     end
% end
% 
% % close all;
% boxColor = repmat('k',4,1);
% 
% figure;
% boxplot(boxData,'symbol','','color',boxColor);
% % set(gca,'ylim',0.45*[0 1],'ytick',[0,0.2,0.4]);
% formatBoxPlot();
% clc
% 
% saveDir = 'C:\Users\Eyal\Google Drive\Research\Publications\14 Neutrophil3D\manuscript\figures\boxAndWhiskers\';
% % export_fig(gcf,[saveDir, 'CI'],'-pdf','-q100','-nocrop');
% % print(gcf,'-dpdf','-painters',[saveDir, 'CI']);

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

minEigE = cell(length(folderDir),1);
maxEigE = cell(length(folderDir),1);

for i = 1:length(folderDir)
    cd(fullfile(path,folderDir{i}));
    disp([num2str(i),'  ', folderDir{i}]);
    
    load('resultsDVC.mat','u','um2vxl','dm');
    load('resultsCell.mat');
    
    nTime   = min([length(u), length(f) - 1]);       % number of DVC increments points
    % -----------------------------------------------------------------------------------------------
    % calculate cell centroid and velocity
    triC       = zeros(nTime + 1,3); % centroid
    
    for j = 1:nTime + 1,
        triC(j,:)           = triCentroid(f{j},v{j});
    end
    dt       = 2;                                                           % minutes
    dCdt{i}        = diff(triC)/dt;
    
    
    % -----------------------------------------------------------------------------------------------
    % -----------------------------------------------------------------------------------------------
    
    nBins = 50;
    bins = linspace(-1,1, nBins + 2);
    bins = bins(2:end-1);
    
    uTri      = cell(nTime,1);       % lagrangian strain interpolated onto surface
    ETri     = cell(nTime,1);
    
    u0 = u;
    Fij = calculateFij(u0,4,'optimal9');
    Eij = calculateEij(Fij);
    
    
    
    for j = 1:nTime
        for k = 1:3, u{j}{k} = u0{j}{k}*um2vxl(k); end
        
        fp = triCenters(f{j},v{j});
        uTri{j}           = triInterp(fp,u{j}(:,1:3),dm);
        
        ETri{j}           = triInterp(fp,Eij{j},dm);
        for ii = 1:size(ETri{j},3)
            eigValETri{j}(ii,:) = eig(ETri{j}(:,:,ii));
        end
        
        m = bsxfun(@minus, fp, mean(fp));
        mNorm = vecNorm(m,2);
        dCdtNorm = dCdt{i}(j,:)/norm(dCdt{i}(j,:));
        
        %   binTri = sum(bsxfun(@times, m , dCdtNorm),2);
        %   bins{i}(:,j) = linspace(min(binTri),max(binTri), nBins)
        binTri = sum(bsxfun(@times, mNorm, dCdtNorm),2);

        
        [~,binInd] = histc(binTri, bins);
        

        
        for ii = 1:nBins
            eigValETri_ = eigValETri{j}(binInd == ii,:);
            
            if isempty(eigValETri_),
                maxEigE{i}(ii,j) = nan;
                minEigE{i}(ii,j) = nan;
            else
                maxEigE{i}(ii,j) = max(eigValETri_(:,3));
                minEigE{i}(ii,j) = min(eigValETri_(:,1));
            end
            
        end
   
    end
    
    
    
    
    %     mechanicsE{i} = calculateMeanMetrics(f,v,uTri);
    
end
%%

minEig = minEigE;
maxEig = maxEigE;

figure

for i = 1:length(folderDir)
    
%     binsNorm{i} = bsxfun(@rdivide, height{i}, height{i}(end,:) - height{i}(1,:));
%     
%     binInt{i} = linspace(max(binsNorm{i}(1,:)),min(binsNorm{i}(end,:)),nBins);
%     
%     for j = 1:size(binsNorm{i},2)
%         
%         maxEig{i}(:,j) = interp1(binsNorm{i}(:,j),maxEPoly{i}(:,j), binInt{i});
%         minEig{i}(:,j) = interp1(binsNorm{i}(:,j),minEPoly{i}(:,j), binInt{i});
%         
%     end
    

    
    plotRange = 1.1*[min(prctile( minEig{i},25,2)) max(prctile(maxEig{i},75,2))];
    
    U1 = inpaint_nans(prctile(maxEig{i},75,2)');
    U2 = inpaint_nans(prctile(maxEig{i},50,2)');
    U3 = inpaint_nans(prctile(maxEig{i},25,2)');
    
    L1 = inpaint_nans(prctile(minEig{i},75,2)');
    L2 = inpaint_nans(prctile(minEig{i},50,2)');
    L3 = inpaint_nans(prctile(minEig{i},25,2)');
    
   
    
    
    
    
    subplot(4,4,i)
    hold on;
    plot(bins, U2, 'r-','linewidth', 3)
    plot(bins, L2, 'b-','linewidth', 3)
    plot([0 0], [-1 1], 'k--', 'linewidth', 2);
    jbfill(bins,U1, U3,'r');
    jbfill(bins,L1, L3,'b');
    hold off
    
    axis([-1,1,plotRange(1) plotRange(2)])
    
    fontSize = 6;
    xlabel('Normalized Distance along velocity axis','FontSize',fontSize);
    ylabel('Principal Strain','FontSize',fontSize);
    title(folderDir{i},'FontSize',fontSize);
    set(gca,'fontsize',fontSize)
    hl = legend('max(Ei)','max(Ei)');
    set(hl,'FontSize',fontSize);
    box on;
    drawnow
    
    
    
end


























% dataInfo.varnames = {'x1','x2','x3','mEmTri','dudtTri'};
% dataInfo.zonename = 'Cell_';
% TPData = convert2Tecplot('FESurfaces', [f(1:end-1),v(1:end-1)],  [mEmTri, dudtTri], dataInfo);
% mat2tecplot(TPData,'mEm_dudt.plt');



%                 figure;
%                 hold on;
%                 plot(binPlanes,maxEPoly,'b.-')
%                 plot(binPlanes,minEPoly,'r.-')
%                 hold off;

%         [vPoly, plane] = calculateContourLines(f{j},v{j},dCdt{i}(j,:),nPoly);
%         EPoly = polyInterp(vPoly, Eij{j},dm);
%         for ii = 1:length(EPoly)
%
%             if all(isnan(EPoly{ii}(:))), EPoly{ii} = 0; end
%
%             eigValEPoly = zeros(size(EPoly{ii},3),3);
%             for jj = 1:size(EPoly{ii},3)
%                 eigValEPoly(jj,:) = eig(EPoly{ii}(:,:,jj));
%             end
%
%             if size(eigValEPoly,1) == 1,
%                 maxEPoly(ii,j) = nan;
%                 minEPoly(ii,j) = nan;
%                 maxEPolyLoc(ii,:) = nan;
%                 minEPolyLoc(ii,:) = nan;
%             else
%
%                 [maxEPoly(ii,j), maxEPolyidx] = max(eigValEPoly(:,3));
%                 [minEPoly(ii,j), minEPolyidx] = min(eigValEPoly(:,1));
%                 maxEPolyLoc(ii,:) = vPoly{ii}(maxEPolyidx,:);
%                 minEPolyLoc(ii,:) = vPoly{ii}(minEPolyidx,:);
%             end
%         end

%         maxEPoly(:,2) = plane(:,4)/(max(plane(:,4)) - min(plane(:,4)));
%         minEPoly(:,2) = plane(:,4)/(max(plane(:,4)) - min(plane(:,4)));
%         zPlane(:,j) = plane(:,4);

%
%         figure;
%         hold on;
%         plot(minEPoly(:,1),minEPoly(:,2),'b.-')
%         plot(maxEPoly(:,1),maxEPoly(:,2),'r.-')
%         hold off;

% close; clear;
% load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsCentroid.mat');
% % http://web.mit.edu/course/other/beh.410/www/Handouts/motility.pdf
% % http://onlinelibrary.wiley.com/doi/10.1002/aic.690391210/pdf
% % http://climate.indiana.edu/RobesonPubs/hanson_etal_vector_correlation.pdf
% 
% histVar = cell(1,length(centroid));
% for i = 1:length(centroid)
%     t = 0:2:2*(size(centroid{i}.CA,1) - 1);
%     histVar{i} = [t' centroid{i}.CA]; 
% end
% 
% histData = cell(1,2);
% tracksT =     histVar(4:6)';   % taxis Naive
% tracksK =     histVar(12:15)'; % Kinesis Naive
% 
% %% T
% nDims = 3;
% msdT = msdanalyzer(nDims, 'µm', 'min');
% msdT = msdT.addAll(tracksT);
% 
% msdT = msdT.computeMSD;
% figure
% msdT.plotMeanMSD(gca, true);
% 
% A = msdT.getMeanMSD;
% % A(isnan(sum(A,2)),:) = [];
% t = A(:, 1); % delay vector
% msd = A(:,2); % msd
% std_msd = A(:,3); % we will use inverse of the std as weights for the fit
% std_msd(1) = std_msd(2); % avoid infinity weight
% 
% ft = fittype('3*a^2*b*(x - b*(1-exp(-x/b)))');
% [fo, gof] = fit(t, msd, ft, 'StartPoint', [1.3 0.4])
% 
% %% K
% 
% nDims = 3;
% msdK = msdanalyzer(nDims, 'µm', 'min');
% msdK = msdK.addAll(tracksK);
% 
% msdK = msdK.computeMSD;
% figure
% msdK.plotMeanMSD(gca, true);
% 
% 
% A = msdK.getMeanMSD;
% % A(isnan(sum(A,2)),:) = [];
% t = A(:, 1); % delay vector
% msd = A(:,2); % msd
% std_msd = A(:,3); % we will use inverse of the std as weights for the fit
% std_msd(1) = std_msd(2); % avoid infinity weight
% 
% ft = fittype('3*a^2*b*(x - b*(1-exp(-x/b)))');
% [fo, gof] = fit(t, msd, ft, 'StartPoint', [1.3 0.4])
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
load('C:\Users\Eyal\Dropbox\Neutrophil3D\data\resultsCentroid.mat');

histVar = cell(1,length(centroid));
for i = 1:length(centroid)
    t = 0:2:2*(size(centroid{i}.CA,1) - 1);
    histVar{i} = [t' centroid{i}.CA]; 
end


histData = cell(1,4);
histData{1} =     histVar(4:6)';   % taxis Naive
histData{2} =     histVar(12:15)'; % Kinesis Naive
histData{3} =     histVar(1:3)';   % taxis LPS
histData{4} =     histVar(7:11)';  % Kinesis LPS

nDims = 3;
figure;
for i = 1:length(histVar)
ma = msdanalyzer(nDims, 'µm', 'min');
ma = ma.addAll(histVar(i));
ma = ma.computeMSD;

A = ma.getMeanMSD;
% A(isnan(sum(A,2)),:) = [];
t = A(:, 1); % delay vector
msd = A(:,2); % msd
std_msd = A(:,3); % we will use inverse of the std as weights for the fit
std_msd(1) = std_msd(2); % avoid infinity weight

options = fitoptions('Method','NonlinearLeastSquares');
options.Lower = [0 0];
options.StartPoint = [1 1];
% options.Display = 'iter';



ft = fittype('3*a^2*b*(x - b*(1-exp(-x/b)))');
[fo, gof] = fit(t, msd, ft, options);

pTime(i,1) = fo.b;
pSpeed(i,1) = fo.a;
rsquare(i,1) = gof.adjrsquare;


subplot(4,4,i); 
msdFit = 3*fo.a^2*fo.b*(t - fo.b*(1-exp(-t/fo.b)));
plot(t,msd,'r^',t,msdFit);

dx = abs(sum(diff(histVar{i}(:,2))));
L = norm(sum(diff(histVar{i}(:,2:4))));

dt = 2;
factor(i,1) = (1 - (dt/pTime(i,1))^(-1)*(1 - exp(-dt/pTime(i,1))));
CI(i,1) = dx/L%*factor(i,1);

end

% pTime(rsquare < 0.95) = nan;
% pSpeed(rsquare < 0.95) = nan;
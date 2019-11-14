histVar = cell(1,length(mechanics));
for i = 1:length(mechanics)
    histVar{i} = (mechanics{i}.F);
end

histData = cell(1,4);
histData{1} =     histVar(12:15)'; % Kinesis Naive
histData{2} =     histVar(4:6)';   % taxis Naive
histData{3} =     histVar(7:11)';  % Kinesis LPS
histData{4} =     histVar(1:3)';   % taxis LPS

Ftotal = cell(1,4);
Jtotal = cell(1,4);
J_ = cell(1,4);

for i=1:4
    for j=1:numel(histData{i})
            Ftotal{i}{j}(:,:) = 2*eye(3)-histData{i}{j}(:,:,1); %flips the normals
            J_{i}{j}(1) = det(2*eye(3)-histData{i}{j}(:,:,1));
        for k = 2:length(histData{i}{j})
            Ftotal{i}{j}(:,:) = (2*eye(3)-histData{i}{j}(:,:,k))*Ftotal{i}{j}(:,:);
            J_{i}{j}(k) = det(2*eye(3)-histData{i}{j}(:,:,k));
        end
        Jtotal{i}{j} = det(Ftotal{i}{j}(:,:));
    end
end

histData = J_;
histData = cellfun(@cell2mat, histData,'uniformoutput',0);

boxData = nan(max(cellfun(@(x) size(x,1), histData)),4);

k = 1;
for i = 1:4
    boxData(1:length(histData{i}),k) = histData{i};
    k = k + 1;
end
boxData(boxData == 0) = NaN;

% boxColor = repmat('k',4,1);
h = figure;
boxplot(boxData,'symbol','', 'widths', 0.5);
set(gca,'ylim',[0.8 1.20],'ytick',[0.80 .9 1 1.1 1.2]);
formatBoxPlot([0.8, 1]);


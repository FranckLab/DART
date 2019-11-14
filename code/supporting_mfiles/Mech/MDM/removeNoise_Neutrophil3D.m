function [data, MEAN, STD, P] = removeNoise_Neutrophil3D(data0,nSTD)
if nargin < 2, nSTD = 2; end

if ~iscell(data0), data0 = {data0}; end


data    = cell(size(data0));
STD     = zeros(size(data0));
MEAN    = zeros(size(data0));
P       = zeros(size(data0));


for i = 1:size(data0,1)
    for j = 1:size(data0,2)
        
        data_       = data0{i,j}(:);
        
        STD(i,j)    = std(data_);
        MEAN(i,j)   = mean(data_);
        [~, P(i,j)] = kstest((data_ - mean(data_))/ std(data_));
        
        noiseFloor = nSTD*STD(i);
        data{i,j}   = data0{i,j} - MEAN(i);
        data{i,j}(  abs(data{i,j}) < abs(noiseFloor)  ) = 0;
        
    end
end



if length(data) == 1, data = data{1}; end



% from Schwann TFM
%     I{i} = I0{i} - mean(noise);
%     noiseLevel{i} = 2*std(noise);
%
%     I{i}(abs(I{i})<abs(noiseLevel{i})) = 0;

end
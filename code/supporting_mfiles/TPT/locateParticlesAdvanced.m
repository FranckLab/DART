%% LocalizeParticles
clc; clear all; close all;

%% User inputs

folder = 'C:\Users\francklab\Desktop\TPTEpiRunFiles\files';
fileformat = '*.mat';
saveformat = 'position_bead%0.3d.mat';

% Check bead identification parameters
thres = 0.5;
minSize = 4;
maxSize = 100;

% Window size
winSize = [7,7,7];

% Figure setting
maxhist = 100;

%% Check user parameters

files = dir(fullfile(folder, fileformat));
load(fullfile(folder, files(2).name))

I = vol{1};
I = I/max(I(:));
BW = I>thres;

CC = bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);
idx = numPixels<3;
idx = find(idx);
for i = 1:length(idx)
    BW(CC.PixelIdxList{idx(i)}) = 0;
end
I = imgaussfilt3(I,0.75);
Im = -I;
Im(~BW) = Inf;
L = watershed(Im);
L(~BW) = 0;

numPixels = regionprops(L, 'Area');
numPixels = double(struct2dataset(numPixels));
figure;hist(numPixels(numPixels<maxhist),30)

%% Find particles in all the images and save the positions. 
beadParam{1}.winSize = winSize;
beadParameter = parseBeadParameter(beadParameter);
for i = 1:length(files)
    
    % Load image and threshold
    load(fullfile(folder, files(2).name))
    I = vol{1};
    I = I/max(I(:));
    BW = I>thres;
    
    % Remove noise and find connect component from watershed
    CC = bwconncomp(BW);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    idx = numPixels<3;
    idx = find(idx);
    for i = 1:length(idx)
        BW(CC.PixelIdxList{idx(i)}) = 0;
    end
    I = imgaussfilt3(I,0.75);
    Im = -I;
    Im(~BW) = Inf;
    L = watershed(Im);
    L(~BW) = 0;
    
    % Create new connected component for confocal images
    numPixels = regionprops(L, 'Area');
    numPixels = double(struct2dataset(numPixels));
    beadBlob = numPixels>minSize & numPixels<maxSize;
    
    % Find centroid of all the eligible blob;
    S = regionprops(L,'Centroid');
    blobPts = round(double(struct2dataset(S)));
    blobPts = blobPts(beadBlob,:);
    
    % Convert to m,n,o coordinates
    blobPts(:,[1,2,3]) = blobPts(:,[2,1,3]);
    x = blobPts;
    
    % Localize particles
    x = radialcenter3dvec(I, x, beadParameter{1});
    
    % Save positions
    file = sprintf(saveformat, i);
    save_name = fullfile(folder, file);
    save(save_name, 'x')
    
    fprintf(strcat('Particles found in file no: %d and saved as:', file, '\n'), i);
end








%%
function beadParameter = parseBeadParameter(beadParameter)

% Define default values
thres = 0.5;
minSize = 5;
maxSize = Inf;
winSize = [7,7,7];
dccd = [1,1,1];
abc = [1,1,1];
forloop = 1;
randNoise = 1/10^7; % Something small
xy = 1; % Assume symmetrical if not given
z = 1;  % Assume symmetrical if not givenf
diameter = 5;   % Assume if not given


for i = 1:length(beadParameter)
    
    p = inputParser;
    addParameter(p,'thres',thres);
    addParameter(p,'minSize',minSize);
    addParameter(p,'maxSize',maxSize);
    addParameter(p,'winSize',winSize);
    addParameter(p,'dccd',dccd);
    addParameter(p,'abc',abc);
    addParameter(p,'forloop',forloop);
    addParameter(p,'randNoise',randNoise);
    addParameter(p,'xy',xy);
    addParameter(p,'z',z);
    addParameter(p,'diameter',diameter);
    
    parse(p,beadParameter{i})
    beadParameter{i} = p.Results;
    
end

end

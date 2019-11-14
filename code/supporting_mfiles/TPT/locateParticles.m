function [x] = locateParticles(I, beadParameter)
% [x] = locateParticles(I, beadParameter) locates particles in the image
%
% INPUTS
% -------------------------------------------------------------------------
%   I:              Input volumetric image
%   beadParameter:  Parameters to detect particle position in images
%
% OUTPUTS
% -------------------------------------------------------------------------
%   x:              Voxel-level estimate of particle center in MxNxO format
%

% Parameters
thres = beadParameter.thres;    %Threshold value
minPixels = beadParameter.minSize;  %Minimum pixel count in blob for bead
maxPixels = beadParameter.maxSize;  %Maximum pixel count in blob for bead

%Process image to help in localizing particles
I = process_images(I);

% Image thresholding
BW = I>thres;

% % % Find bead blobs
% % CC = bwconncomp(BW);
% % numPixels = cellfun(@numel,CC.PixelIdxList);
% % beadBlob = numPixels>minPixels & numPixels<maxPixels;
% % 
% % % Find centroid of all the eligible blob;
% % S = regionprops(CC,'Centroid');
% % blobPts = round(double(struct2dataset(S)));
% % blobPts = blobPts(beadBlob,:);
% % temp = blobPts;
% % toc

% Remove noise and find connect component from watershed
CC = bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);
idx = numPixels<minPixels;
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
beadBlob = numPixels>minPixels & numPixels<maxPixels;

% Find centroid of all the eligible blob;
S = regionprops(L,'Centroid');
blobPts = round(double(struct2dataset(S)));
blobPts = blobPts(beadBlob,:);


% Convert to m,n,o coordinates
blobPts(:,[1,2,3]) = blobPts(:,[2,1,3]);
x = blobPts;

end


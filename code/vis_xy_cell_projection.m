%% Visualize xy cell projection

clear; close all;

% Visualization data file
file = '20190101_A01_dispvis01.mat';
load(file)

% Convert label to rgb after xy projection
BW = max(BW,[],3);
I = label2rgb(BW,'lines', 'k');

% Add cell number to the I according to position
for i = 1:max(BW)
    stats = regionprops(BW == i, 'Centroid');
    cell_centroid = stats.Centroid(1:2);
    
    I = insertText(I, cell_centroid, num2str(i), ...
        'FontSize', 72, 'BoxOpacity', 0, 'TextColor', 'white',...
        'AnchorPoint', 'LeftTop', 'Font', 'Arial');
end

% Visualize
figure;
imshow(I)


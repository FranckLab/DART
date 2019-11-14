function [] = run_filter_cell_img(i)


addpath(genpath('./MatlabCode'))

% Folder to find cell surface in
folder = './20181028';

% Find valid MP dir in the folder
not_valid_dir = {'.','..','SDS'};
subdir = dir(folder);
subdir = subdir([subdir.isdir]);
valid_dir = false(1, length(subdir));
for j = 1:length(valid_dir)    
    if sum(strcmp(subdir(j).name, not_valid_dir)) == 0
        valid_dir(j) = 1;
    end
end
subdir = subdir(valid_dir & [subdir.isdir]);
subdir = {subdir.name};

% File format for cell image file
cell_file = '*cell*.mat';
files = dir(fullfile(folder, subdir{i}, cell_file));

% Run filter_cellimg for individual cell image files
for j = 2:length(files)
    filename = fullfile(files(j).folder, files(j).name);
    filter_cell_img(filename)
end

end
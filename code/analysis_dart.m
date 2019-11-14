function [] = analysis_dart(folder, i)
% Code to run cell segmentation and dart analysis

tic
addpath(genpath('./MatlabCode'))

% Find the valid MP dir in the folder
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

% Create table to store data. Note that we are computing many mechanical
% quantities which have not been used in this paper. 
varNames = {'date', 'well', 'exp_start_time', 'induce_type', 'drug', 'conc', ...
    'cell', 'img_time', 'is_disp_checked', 'volume', 'area', ...
    'J', 'eigValU', 'axisAngleR', 'max_u', 'avg_u', 'cell_center', ...
    'shape_anisotropy', 'perimeter_2d', 'projected_2dtarea', 'in_avg_u',...
    'in_med_u', 'D', 't1', 't2', 't3', 't4', 't5', 'c1', 'c2', 'c3', 'c4', 'c5', ...
    'disp_vol', 'pull_', 'pull_rms', 'pull_peak', 'push_', 'push_rms', 'push_peak',...
    'shear_', 'shear_rms', 'shear_peak', 'equiv_dia', 'max_ferret_dia'};
mech_table = cell2table(cell(0,length(varNames)), 'VariableNames', varNames);

% Load condition file
data_condition = importConditionFile(fullfile(folder, 'Conditions.xlsx'));
data_row = cell(1, length(varNames));
[data_row{1:6}] = data_condition{i,:};

% T-PT displacement results
tpt_file = strcat('resultsTPT_', subdir{i}, '.mat');
u = load(fullfile(folder, subdir{i}, tpt_file));
x = u.x; tracks = u.track; clear u;

% Threshold values for cell cluster segmentation. These values need to be
% manually determined. 
thres = importfile(fullfile(folder, 'thres.xlsx'));

for t = 2:2%Select the timepoint for which the analysis should be run
    
    % Prepare the t-pt displacement
    u.x0 = x{1}{1};
    u.x1 = x{t+1}{1};
    u.track = tracks{t}{1};
    calib = [0.32,0.32,0.6];
    
    % Load and segement cell image
    vol_file = fullfile(files(t).folder, files(t).name);
    [fv, label_BW] = compute_cell_surface(vol_file, thres.thres(i), calib);
                            
    display(sprintf('Total number of cells: %d', length(fv)));
    toc
    for j = 1:length(fv)
        
        display(sprintf('Current cell: %d', j));
        
        % Fill data_row;
        data_row{7} = j; %cell
        data_row{8} = t; %time
        data_row(9:end) = {NaN};
        
        % Compute MDM quantities
        [mdm] = compute_mech(fv{j}, u, calib, size(label_BW));
        
        % Compute other mech quantities
        [dart, pts, u_pts] = compute_dart(label_BW,j, u, calib);
        data = [mdm, dart];
        [data_row{10:end}] = data{:};
        mech_table = [mech_table; data_row];
        
        % Save data for visualization purposes
        vis{j} = {pts, u_pts, fv{j}, data_row{17}, calib};
        
        toc
    end
    
    % Save data for visualization purposes
    save_file = strcat(folder(3:end), '_', subdir{i}, '_', ...
                sprintf('dispvis%0.2d', t), '.mat');
    save_folder = files(1).folder;
    save(fullfile(save_folder, save_file), 'vis', 'label_BW');
end

save_table = fullfile(folder,  subdir{i}, 'dart.csv');
writetable(mech_table, save_table);

display('Analysis done! :)')

end

function [fv, label_BW] = compute_cell_surface(vol_file, thres, calib)

% file load image
load(vol_file, 'vol')
vol = double(vol);
vol = vol/(2^12);
I = vol;

%% Uncomment this section and only run it to figure out the thres value. 
% Basically, set the code to stop at the end of this section, and run this
% section with multiple thres value. Select the thres value which segements
% the cell cluster the best. Later comment this section back. Also, try to 
% be consistent in your criteria for evaluating what is a good thres value.
% This step takes a long time. Wish you luck and patience. 
 
% % % thres = 0.45; %trial thres value
% % % 
% % % % Segment the cell cluster
% % % vol = imbinarize(I, adaptthresh(I, thres));
% % % minPixels = 5000;
% % % vol = bwareaopen(vol, minPixels, 6);
% % % % vol = bwdist_aspect(vol, calib); % Too computational expensive. Computes
% % % % bwdist taking into account the different aspect along x, y, and z axis
% % % vol = bwdist(vol); % So will run this version instead, less correct but faster
% % % vol = vol<3;
% % % minPixels = 8000; % Remove small image volumes
% % % vol = bwareaopen(vol, 8000, 6);
% % % vol = imfill(vol,'holes'); % Computational expensive so don't run it over here
% % % 
% % % % Visualize to evaluate segementation performance
% % % figure; 
% % % imshow3D(cat(2,I*6,vol)); 

%% Once the threshold value has been found, comment the section in the top, and run this instead

vol = imbinarize(I, adaptthresh(I, thres));
minPixels = 5000; % Remove small image volumes
vol = bwareaopen(vol, minPixels, 6);
vol = bwdist_aspect(vol, calib);
vol = vol<3;
minPixels = 8000; % Remove small image volumes
vol = bwareaopen(vol, 8000, 6);
vol = imfill(vol,'holes');

%% Further processing 

% Remove cell region segmented from top and bottom of the image
vol(:,:,1:8) = 0;
vol(:,:,114:end) = 0;

% Create label matrix
label_BW = labelmatrix(bwconncomp(vol)); clear vol

% Find cell surface
for i = 1:max(label_BW(:))
    
    %Segement the selected cell
    BW = label_BW == i;
    
    % Grow cell by dist of 5;
    BW = bwdist_aspect(BW, calib);
    BW = BW < 5;
    
    % Crop top and bottom cell segmentation
    BW(:,:,1:8) = 0;
    BW(:,:,114:end) = 0;
    fv{i} = smoothpatch(isosurface(BW, 0.5), 1, 50);
    fv{i}.vertices = fv{i}.vertices(:,[2,1,3]); % fix matlab x,y flip
end

end


function [data, pts, u_pts] = compute_dart(label_BW, cell_idx,  u, calib)
data={};

% Prepare displacemnt
x0 = u.x0(u.track>0, :);
x1 = u.x1(u.track(u.track>0), :);
u = x1-x0;
x = num2cell(x1, 1);
u = removeRigidDriftRotation(u, x);
u = u.*calib;
x = num2cell(x1.*calib, 1);
x1 = x1.*calib;

% Create cell binary masks
BW = label_BW == cell_idx;
other_cell_BW = label_BW>0; other_cell_BW(BW) = 0;

% New addition
cc = regionprops3(BW, 'centroid');
center = cc.Centroid;


% Uniformly sample the image space. This is the spacing at which we will
% compute the displacement
dm = 12;
for i = 1:3
    m{i} = 1:dm:size(BW, i);
end
[m{:}] = ndgrid(m{:});
idx_ = sub2ind(size(BW), m{1}(:), m{2}(:), m{3}(:));
u_BW = zeros(size(BW));
u_BW(idx_) = 1;

% Consider valid points around the cell
BW = bwdist_aspect(BW, calib);
BW = BW < 70 & BW>0;
BW(:,:,1:8) = 0;
BW(:,:,114:end) = 0;
BW(other_cell_BW)=0; % Remove BW covering outher cells

% Points at which displacement will be computed
IDX = find(BW.*u_BW);
[pts(:,1), pts(:,2), pts(:,3)] = ind2sub(size(BW), IDX);
pts = pts.*calib;

% Smoothen and interpolate displacement at the given points
[idx, dist] = knnsearch(x1, x1, 'k', 5);
mu = 0;
sigma = 3;
weight = pdf('Normal', dist, mu, sigma);
weight = weight./sum(weight,2);
idx_shape = size(idx);
for j = 1:3
    u(:,j) = sum(reshape(u(idx, j), idx_shape).*weight, 2);
    F = scatteredInterpolant(x{:}, u(:,j), 'linear', 'none');
    u_pts(:,j) = F(pts);
end

% remove nan points
idx = all(~isnan(u_pts),2);

%%% Displacement is now ready for further analysis
temp_results = {};


% Find cell center
BW = label_BW == cell_idx;
CC = bwconncomp(BW);

% Cell 3D assymetry
s = regionprops3(CC, 'Centroid', 'EigenValues');
cell_center = s.Centroid([2,1,3]).*calib;
data{end+1} = cell_center;

% Cell Project 2D Assymetry
CC = bwconncomp(sum(BW,3)>0);
s = regionprops(CC, 'Eccentricity', 'MajorAxisLength', 'MinorAxisLength',...
    'Perimeter', 'FilledArea', 'EquivDiameter', 'MaxFeretProperties');
data{end+1} = s.MajorAxisLength/s.MinorAxisLength;
data{end+1} = s.Perimeter*calib(1);
data{end+1} = s.FilledArea*prod(calib(1:2));
equiv_dia = s.EquivDiameter;
max_feret_dia = s.MaxFeretDiameter;

% Average inward displacement magnitude
r = pts - cell_center;
norm_vec = r./sqrt(sum(r.^2,2));
in_u = sum(u_pts.*norm_vec, 2);
shear_u = sqrt(sum(u_pts.^2, 2)) - abs(in_u);
data{end+1} = nanmean(in_u);
data{end+1} = nanmedian(in_u);

% Second moment of inward displacement
dv = 12^3*prod(calib);
M_uin = zeros(3,3);
for i =1:3
    for j = 1:3
        M_uin(i,j) = nansum(dv*in_u.*r(:,i).*r(:,j));
        M_uin(i,j) = nansum(dv*in_u.*norm_vec(:,i).*norm_vec(:,j));
    end
end
[~,D] = eig(M_uin);
data{end+1} = [diag(D)]';

% DART computation 
uin_BW = zeros(size(BW));
uin_BW(IDX) = in_u;
uin_BW_reduced = zeros(size(m{1}));
uin_BW_reduced(:) = uin_BW(idx_);

ushear_BW = zeros(size(BW));
ushear_BW(IDX) = shear_u;
ushear_BW_reduced = zeros(size(m{1}));
ushear_BW_reduced(:) = ushear_BW(idx_);

disp_thres = 0.4;

CC = bwconncomp(uin_BW_reduced>disp_thres);
s = regionprops3(CC, uin_BW_reduced, 'Volume', 'MaxIntensity', 'MeanIntensity');
s = sortrows(s, 'Volume', 'descend');
s = [[s.Volume], [s.MaxIntensity], [s.MeanIntensity]];
for i = 1:5
    if i < size(s, 2)
        if s(1, i)<11
            data{end+1} = NaN(1,3);
        else
            data{end+1} = s(i,:);
        end
    else
        data{end+1} = NaN(1,3);
    end
end

CC = bwconncomp(uin_BW_reduced<-disp_thres);
s = regionprops3(CC, -uin_BW_reduced, 'Volume', 'MaxIntensity', 'MeanIntensity');
s = sortrows(s, 'Volume', 'descend');
s = [[s.Volume], -[s.MaxIntensity], -[s.MeanIntensity]];
for i = 1:5
    if i < size(s, 2)
        if s(1, i)<11
            data{end+1} = NaN(1,3);
        else
            data{end+1} = s(i,:);
        end
    else
        data{end+1} = NaN(1,3);
    end
end

data{end+1} = sum(~isnan(in_u));

% push, pull and shear displacements
push = bwareaopen(uin_BW_reduced>disp_thres, 12);
pull = bwareaopen(uin_BW_reduced<-disp_thres, 12);
shear = bwareaopen(ushear_BW_reduced>disp_thres, 12);

for i = 1:3
    m{i} = 1:dm:size(BW, i);
end
[m{:}] = meshgrid(m{:});
for i = 1:3
    m{i} = m{i} - center(i);
end
[THETA,RHO,R] = cart2sph(m{1},m{2},m{3});

push_ = [];
pull_ = [];
push_rms = [];
pull_rms = [];
push_peak = [];
pull_peak = [];
shear_ = [];
shear_rms = [];
shear_peak = [];

for i = 1:2
    if i == 1
        IDX_RHO = RHO>=0;
    end
    if i == 2
        IDX_RHO = RHO<0;
    end
    for j = 1:8
        theta_min = -pi+(j-1)*pi/4;
        theta_max = -pi+j*pi/4;
        IDX = IDX_RHO & THETA>theta_min & THETA <= theta_max;
        
        pull_(end+1) = sum(pull(IDX));
        push_(end+1) = sum(push(IDX));
        shear_(end+1) = sum(shear(IDX));
        
        if pull_(end)>0
            temp = -uin_BW_reduced(IDX & pull);
            pull_rms(end+1) = rms(temp(temp>0));
            pull_peak(end+1) = max(temp(temp>0));
        else
            pull_rms(end+1) = NaN;
            pull_peak(end+1) = NaN;
        end
        
        if push_(end)>0
            temp = uin_BW_reduced(IDX & push);
            push_rms(end+1) = rms(temp(temp>0));
            push_peak(end+1) = max(temp(temp>0));
        else
            push_rms(end+1) = NaN;
            push_peak(end+1) = NaN;
        end
        
        if shear_(end)>0
            temp = ushear_BW_reduced(IDX & shear);
            shear_rms(end+1) = rms(temp(temp>0));
            shear_peak(end+1) = max(temp(temp>0));
        else
            shear_rms(end+1) = NaN;
            shear_peak(end+1) = NaN;
        end
        
    end
end

data{end+1} = pull_;
data{end+1} = pull_rms;
data{end+1} = pull_peak;

data{end+1} = push_;
data{end+1} = push_rms;
data{end+1} = push_peak;

data{end+1} = shear_;
data{end+1} = shear_rms;
data{end+1} = shear_peak;

data{end+1} = equiv_dia;
data{end+1} = max_feret_dia;

end


function D = bwdist_aspect(BW, calib)

% Create original grid
size_calib = size(BW).*calib - calib;
for i = 1:3
    n{i} = 0:calib(i):size_calib(i);
end
[n{:}] = ndgrid(n{:});

% Create equal grid
new_calib = [calib(1),calib(1),calib(1)];
for i = 1:3
    m{i} = 0:new_calib(i):size_calib(i)+new_calib(i);
end
[m{:}] = ndgrid(m{:});

% Interpolate with nearest neighbor to equal grid
BW = interpn(n{:}, single(BW), m{:}, 'nearest');

% Computer BW dist
D = bwdist(BW>0);

% Convert D to original grid
D = interpn(m{:}, D, n{:}, 'linear');
end



function [data] = plot_coneplot(x, u, fv, cell_center, calib)

% Calibrate fv
fv.vertices = fv.vertices.*calib;
fv.vertices = fv.vertices(:,[2,1,3]); % convert to xyz format

% Find center of the cell
% cell_center = regionprops(BW, 'centroid');
% cell_center = cell_center.Centroid;
% cell_center = cell_center([2,1,3]).*calib

% pt_vect
norm_vec = x - cell_center;
norm_vec = norm_vec./(sqrt(sum(norm_vec.^2, 2)));

% Magnitude aligned along the inward vector
umag = sum(u.*norm_vec,2);
cmax = max(abs(umag));



% Plot figure
fv = smoothpatch(isosurface(BW, 0.5), 1, 50);
fv.vertices = fv.vertices.*calib;

% umag = sqrt(sum(u.^2, 2));
figure;
hold all
hc = coneplot(x(:,2),x(:,1),x(:,3),u(:,2),u(:,1),u(:,3),0.05,'nointerp');
fvc = repmat(umag(:)',[42 1]);
set(hc, 'facecolor', 'flat', 'facevertexcdata', fvc(:))
hc.EdgeColor = 'none';
hc.AmbientStrength = 0.6;
hc.DiffuseStrength = 0.75;
hc.SpecularStrength = 0.4;
hc.SpecularExponent = 3;
% colormap(hot)
colormap(flipud(redblue(100)))
% colorbar
caxis([-cmax cmax])
freezeColors

hs = patch(fv);
hs.FaceColor = [1,1,1]*225/255;
hs.EdgeColor = 'none';
hs.AmbientStrength = 0.40;
hs.DiffuseStrength = 0.50;
hs.SpecularStrength = 0.4;
hs.SpecularExponent = 3;
axis image;
hl = light;
lightangle(hl,160, 20)
view([132.61 6.5]);
camva('manual')
lighting gouraud

h = gcf;
set(gca,'color','none','Visible','off')
set(h,'color',[1,1,1]);

% export_fig(h, 'preinduced2.jpg', '-jpeg', '-opengl')
close all;
data = 0;
end

function tableout = importfile(workbookFile,sheetName,startRow,endRow)
%% Input handling

% If no sheet is specified, read first sheet
if nargin == 1 || isempty(sheetName)
    sheetName = 1;
end

% If row start and end points are not specified, define defaults
if nargin <= 3
    startRow = 1;
    endRow = 48;
end

%% Import the data
[~, ~, raw] = xlsread(workbookFile, sheetName, sprintf('A%d:B%d',startRow(1),endRow(1)));
for block=2:length(startRow)
    [~, ~, tmpRawBlock] = xlsread(workbookFile, sheetName, sprintf('A%d:B%d',startRow(block),endRow(block)));
    raw = [raw;tmpRawBlock]; %#ok<AGROW>
end
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
stringVectors = string(raw(:,1));
stringVectors(ismissing(stringVectors)) = '';
raw = raw(:,2);

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
I = cellfun(@(x) ischar(x), raw);
raw(I) = {NaN};
data = reshape([raw{:}],size(raw));

%% Create table
tableout = table;

%% Allocate imported array to column variable names
tableout.well = stringVectors(:,1);
tableout.thres = data(:,1);
end

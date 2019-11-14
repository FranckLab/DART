function [] = runTPT(i)
% Here, i denothes the ith multi-point sub folder on which T-PT will run
% inside the prescribed folder. 


% Run TPT for multipoints in the given folder
folder = './20181028';

% Add required matlab code
addpath(genpath('./MatlabCode'))

% If running on Brown's ccv, uncomment the following 4 lines of code
% pc = parcluster('local')
% mkdir(strcat('./scratch/', getenv('SLURM_JOB_ID')))
% pc.JobStorageLocation = strcat('./scratch/', getenv('SLURM_JOB_ID'))
% parpool(pc, 6)

% Find valid multipoint dir in the folder
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

%Change path to the subdir
cd(fullfile(folder, subdir{i}))

% TPT_parameters
fileInfo{1}{1} = '*bead*.mat'; %filename
beadParam{1}.thres = 0.075; %0.075
beadParam{1}.minSize = 25;
beadParam{1}.maxSize = 230;
beadParam{1}.winSize = [5, 5, 7];
tptParam{1}.knnFM = 5;
tptParam{1}.fmThres = 2;
tptParam{1}.outlrThres = 5;

% Run TPT
[x, track] = funTPT(fileInfo, beadParam, tptParam);

% Save results
savename = strcat('resultsTPT_', subdir{i})
save(savename, 'x', 'track')
end

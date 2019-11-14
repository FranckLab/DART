function [] = run_nd2mat(i)


%% Preparation

% Addbioformats to path
addpath(genpath('./MatlabCode'))

% folder and file names
folder{1} = './20190101';
folder{2} = './20190101';

filename{1} = '3d tfm 10-28-18 wk1.nd2'; % Images before cell lyses
filename{2} = '3d tfm 10-28-18 wk1 after SDS.nd2'; % Image after cell cell lyses

save_prefix{1} = '20181028';
save_prefix{2} = '20181028_SDS';

% Convert .nd2 to .mat files for ith folder and file
nd2mat(filename{i}, folder{i}, save_prefix{i})

end


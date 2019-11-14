function [] = vis_wholeimg(cell_file, tpt_file, t, ...
                            fig_file, data_file)
%% Visualize displacement field in the whole image around each cell.
% % clear;
seed_frac = 0.00015;
% cell_file = '20181028_E03_cell_T001.mat';
% tpt_file = 'resultsTPT_E03.mat';
% t = 1;


%Cell image
load(cell_file, 'vol');
vol = double(vol);
vol = vol/max(vol(:));
pd = fitdist(vol(:), 'normal');
idx = vol(:) < (pd.mu + 4*pd.std);
pd = fitdist(vol(idx), 'normal');
thres = pd.mu + 3.5*pd.std;

% figure;hist(vol(:), 100);
% Threshold image
vol = vol>thres;
minPixels = 1500;
vol = bwareaopen(vol, minPixels, 6);
% figure; imshow3D(cat(2,I,vol));
%%
vol = imfill(vol, 'holes');
fv = smoothpatch(isosurface(vol, 0.5), 1, 80);
% fv.vertices = fv.vertices(:,[2,1,3]);

%%
% Load TPT
load(tpt_file)
clear u
u.x0 = x{1}{1};
u.x1 = x{t+1}{1};
u.track = track{t}{1};
calib = [0.32,0.32,0.6];
x0 = u.x0(u.track>0, :);
x1 = u.x1(u.track(u.track>0), :);
u = x1-x0;
x = num2cell(x1, 1);
% Remove rigid motion and calibrate displacement
u = removeRigidDriftRotation(u, x);
% Calibrate into proper units
idx = x0(:,3) > 115;
u(idx,:) = 0;
u = u.*calib;
x = x1;


% Compute points around cell to plot cone plots
D = bwdist(vol);
idx = D>1 & D<100;
D = zeros(size(D));
D(idx) = 1;
clear pts
[pts(:,1), pts(:,2), pts(:,3)] = ind2sub(size(D), find(idx));
idx = rand(size(pts,1), 1);
idx = idx<seed_frac;
pts = pts(idx,:);
y = zeros(size(pts));
for i = 1:3
    F = scatteredInterpolant(x(:,1), x(:,2), x(:,3), u(:,i), 'linear', 'none');
    y(:,i) = F(pts);
end
u = y;
x = pts;
umag = sqrt(sum(u.^2,2));
% idx = umag>0.25;
% u = u(idx, :);
% x = x(idx,:);
% umag = umag(idx);



figure;
hold all
hc = coneplot(x(:,2),x(:,1),x(:,3),u(:,2),u(:,1),u(:,3),0.03,'nointerp');
fvc = repmat(umag(:)',[42 1]);
set(hc, 'facecolor', 'flat', 'facevertexcdata', fvc(:))
hc.EdgeColor = 'none';
hc.AmbientStrength = 0.6;
hc.DiffuseStrength = 0.75;
hc.SpecularStrength = 0.4;
hc.SpecularExponent = 3;
colormap(hot)
freezeColors

hs = patch(fv);
hs.FaceColor = [60 60 60]/255;
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
% set(gca,'color',0.9*[1,1,1],'Visible','off')
set(gca,'Visible','off')
set(h,'color','k');
c = colorbar;
c.Color = 0.7*[1,1,1];

savefig(fig_file)
close all;
save(data_file, 'x', 'u', 'fv')

end


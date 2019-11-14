%% Visualize radial displacement around the cell.

% Visualization data file
file = '20190101_A01_dispvis01.mat';
load(file)

% Visualize displacement for 'ith' cell cluster in the well
i = 1;


% Extract data
x = vis{i}{1};
u = vis{i}{2};
fv = vis{i}{3};
cell_center = vis{i}{4};
calib = vis{i}{5};

% Calibrate fv
fv.vertices = fv.vertices.*calib;
fv.vertices = fv.vertices(:,[2,1,3]); % convert to xyz format

% pt_vect
norm_vec = x - cell_center;
norm_vec = norm_vec./(sqrt(sum(norm_vec.^2, 2)));

% Magnitude aligned along the inward vector
umag = sum(u.*norm_vec,2);
cmax = max(abs(umag));


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
colormap(flipud(redblue(100)))
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

%%% Apply this setting before saving the figure
% h = gcf;
% set(h,'color',[1,1,1]);
% colorbar;

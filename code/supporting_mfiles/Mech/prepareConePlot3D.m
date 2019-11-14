function [FV, disp, cmax, cell_center] = prepareConePlot3D(cell_idx, L, ...
                                                        u, calib, cone_vol)


% Padarray in case the cell is on the border
sizeL = size(L);
w = 126;
L = padarray(L, [w, w, 0]); % Watershed segmented labels

% x-y cell center
xy = regionprops(L, 'centroid');
xy = round(xy(cell_idx).Centroid); xy([1,2,3]) = xy([2,1,3]);
cell_center = xy.*calib;
xy(3) = 0;

% Move the particle displacement to the new padded image grid
x0 = u.x0(u.track>0, :);
x1 = u.x1(u.track(u.track>0), :);
u = x1-x0;

% Remove broder points
minxyz = [20,20,12];
maxxyz = sizeL - [20,20,12];
idx = any(x1 < minxyz, 2) | any(x1 > maxxyz, 2);

% Bring the xy grid of tpt displacements points to the grid of the image
x = x1 + [w, w, 0]; clear x0 x1
x = num2cell(x, 1);

% Remove rigid motion and calibrate displacement
u = removeRigidDriftRotation(u, x);
u(idx,:) = 0.00001;
u = u.*calib;

% Crop the image to generate cell surface
BW = double(L(xy(1)-w:xy(1)+w-1, xy(2)-w:xy(2)+w-1, :));
BW_maincell = BW == cell_idx;
BW_othercell = BW.*double(~BW_maincell);
BW_othercell = BW_othercell > 0;
BW = double(~(BW>0));

% Surface for the main cell
BW_maincell = bwareaopen(BW_maincell, 1500, 6);
FV_main = smoothpatch(isosurface(BW_maincell, 0.5), 0, 20,1);
FV_main = reducepatch(FV_main, 2500, 'fast');
FV_main.vertices = FV_main.vertices(:, [2,1,3]); % Convert from xyz to mno format
FV_main.vertices = FV_main.vertices - 2*[w,w,0] + xy;

% Save for ploting
FV = cell(4,1);
FV{1,1} = FV_main.faces-1;
FV{2,1} = FV_main.vertices.*calib;

% Find valid points to plot conepts
[cpts(:,1), cpts(:,2), cpts(:,3)] = ind2sub(size(BW), find(cone_vol.*BW));
% cpts = removeBorderPts(cpts, size(BW), 0);
if size(cpts, 1)>7000
     idx = randperm(size(cpts, 1),7000);
     cpts = cpts(idx, :);
end
conepts_interp = cpts + [xy(1)-w, xy(2)-w, 0];

% Interpolate displacement at cone points
for j = 1:3
    F = scatteredInterpolant(x{:}, u(:,j), 'linear', 'none');
    u_cone(:,j) = F(conepts_interp(:,[1,2,3]));
end
u_cone(isnan(u_cone)) = 0;

% Save displacements
cmag = sqrt(sum(u_cone.^2, 2));
cmax = max(cmag);
disp{1,1} = u_cone;
disp{2,1} = (conepts_interp-[w,w,0]).*calib;

end

function pts = removeBorderPts(pts, sizeI, w)

minxyz = [w, w, 0] + [0, 0, 12];
maxxyz = sizeI - [w, w, 0] - [15, 15, 12];

idx_low = pts < minxyz;
idx_low = sum(idx_low, 2) > 0;

idx_high = pts > maxxyz;
idx_high = sum(idx_high, 2) > 0;

idx = (idx_low + idx_high)>0;
pts = pts(~idx, :);

end

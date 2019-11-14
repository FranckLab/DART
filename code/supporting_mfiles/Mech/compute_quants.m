function [FV, disp, cmax] = compute_quants(cell_idx, L, u, calib, cone_vol)

%%%%% Find cell surface
% Compute centroid
cell_center = regionprops(L, 'centroid');
cell_center = round(cell_center(cell_idx).Centroid); 
cell_center([1,2,3]) = cell_center([2,1,3]);

% Crop the image to generate cell surface
BW_maincell = L == cell_idx;
BW_othercell = double(L).*double(~BW_maincell);
BW_othercell = BW_othercell > 0;

% Surface for the main cell
FV = smoothpatch(isosurface(BW_maincell, 0.5), 1, 5);
% FV_main = reducepatch(FV_main, 2500, 'fast');

% If the cell mask is within valid range
isvalid = 0;
if (min(FV.vertices(:,1)) < 10 | max(FV.vertices(:,1)) > size(L,1)-10 ...
    | min(FV.vertices(:,2)) < 10 | max(FV.vertices(:,2)) > size(L,2)-10 ...
    | min(FV.vertices(:,3)) < 5 | max(FV.vertices(:,3)) > size(L,3)-15)

    isvalid = 1;
end

% Output to 
FV.vertices = FV.vertices.*calib;
fp = triCenters(FV.vertices, FV.vertices);


%%%%% Interpolate displacement to cell surface
% Move the particle displacement to the new padded image grid
x0 = u.x0(u.track>0, :);
x1 = u.x1(u.track(u.track>0), :);
u = x1-x0;
x = num2cell(x1, 1);

% Remove rigid motion and calibrate displacement
u = removeRigidDriftRotation(u, x);
u = u.*calib;

% Interpolate
for j = 1:3
    F = scatteredInterpolant(x{:}, u(:,j), 'linear');
    uTri(:,j) = F(fp);
end

mechanics = computeMeanMetrics(f,v,uTri)
mechanics{end+1} = isvalid;


end


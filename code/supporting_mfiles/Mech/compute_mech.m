function [mechanics] = compute_mech(FV, u, calib, sizeL)


% If the cell mask is within valid range
isvalid = 1;
if (min(FV.vertices(:,1)) < 10 | max(FV.vertices(:,1)) > sizeL(1)-10 ...
    | min(FV.vertices(:,2)) < 10 | max(FV.vertices(:,2)) > sizeL(2)-10 ...
    | min(FV.vertices(:,3)) < 5 | max(FV.vertices(:,3)) > sizeL(3)-15)

    isvalid = 0;
end

% Output to 
FV.vertices = FV.vertices.*calib;
% FV.faces = FV.faces(:,[1,3,2]);
fp = triCenters(FV.faces, FV.vertices);


%%%%% Interpolate displacement to cell surface
% Move the particle displacement to the new padded image grid
x0 = u.x0(u.track>0, :);
x1 = u.x1(u.track(u.track>0), :);
% [x] = remove_rigid_motion(x0, x1, sizeL);
u = x1-x0;

x = num2cell(x1, 1);

% Remove rigid motion and calibrate displacement
u = removeRigidDriftRotation(u, x);

% % Remove broder points
minxyz = [20,20,5];
maxxyz = sizeL - [20,20,10];
idx = any(x1 < minxyz, 2) | any(x1 > maxxyz, 2);
u(idx,:) = 0;

u = u.*calib;
x = num2cell(x1.*calib, 1);
x1 = x1.*calib;

% Smoothen displacement
[idx, dist] = knnsearch(x1, x1, 'k', 5);
mu = 0;
sigma = 2;
weight = pdf('Normal', dist, mu, sigma);
weight = weight./sum(weight,2);
idx_shape = size(idx);


% Interpolate displacement at fp.
for j = 1:3
    u(:,j) = sum(reshape(u(idx, j), idx_shape).*weight, 2);
    F = scatteredInterpolant(x{:}, u(:,j), 'linear', 'none');
    uTri(:,j) = F(fp);
end
uTri(isnan(uTri)) = 0;

% Compute mean metrics
mechanics = computeMeanMetrics(FV.faces,FV.vertices,uTri);
end


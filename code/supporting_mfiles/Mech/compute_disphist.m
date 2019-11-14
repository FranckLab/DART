function [data, temp] = compute_disphist(fv, L, cell_idx, u, calib)
%Compute displacement histogram from the cell surface.


% Parameters
min_volume = 100;
density = 0.0125;


% Load displacement
x0 = u.x0(u.track>0, :);
x1 = u.x1(u.track(u.track>0), :);
u = x1-x0;
x = num2cell(x1, 1);
% Remove rigid motion and calibrate displacement
u = removeRigidDriftRotation(u, x);

% Calibrate into proper units
fv.vertices = fv.vertices.*calib;
% fv.faces = fv.faces(:,[1,3,2]);% flip faces
u = u.*calib;
x = num2cell(x1.*calib, 1);
x1 = x1.*calib;

% Compute cell volume and area
cell_volume = triVolume(fv.faces, fv.vertices);
cell_area = triArea(fv.faces, fv.vertices);


data=NaN;

% Generate random points around the cell in the valid range
min_xyz = [20,20,5].*calib;
max_xyz = (size(L) - [20,20,10]).*calib;
[pts, d, n_vec] = generate_rand_pts(min_xyz, max_xyz, density, fv, L==cell_idx);

% Smoothen displacement
[idx, dist] = knnsearch(x1, x1, 'k', 5);
mu = 0;
sigma = 3.6;
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
pts = pts(idx, :);
d = d(idx, :);
n_vec = n_vec(idx, :);
u_pts = u_pts(idx, :);

% Script to plot cone plot
temp = plot_coneplot(pts, u_pts, fv);

% Compute histogram
bins = 0:1.5:30;

try
    data = compute_bins(bins, d, u_pts, n_vec);
end

end

function data = compute_bins(bins, d, u_pts, n_vec)

% Compute norml and tanganetial displacement vector
u_norm = n_vec.*(dot(u_pts, n_vec, 2));
u_tang = u_pts - u_norm;

% Compute magnitudes
u_mag = sqrt(sum(u_pts.^2 ,2));
un_mag = sqrt(sum(u_norm.^2 ,2));
un = dot(u_pts, n_vec, 2);
ut_mag = sqrt(sum(u_tang.^2 ,2));


% Compute histogram
umag = zeros(length(bins)-1, 1);
unorm = zeros(length(bins)-1, 1);
unorm_vec = zeros(length(bins)-1, 1);
utang = zeros(length(bins)-1, 1);
umag_avg = zeros(length(bins)-1, 1);
unorm_avg = zeros(length(bins)-1, 1);
unorm_vec_avg = zeros(length(bins)-1, 1);
utang_avg = zeros(length(bins)-1, 1);

dist = zeros(length(bins)-1, 1);

for i = 1:length(bins)-1
    dist(i) = mean(bins(i:i+1));
    idx = d>=bins(i) & d<bins(i+1);
    
    temp = maxk(u_mag(idx), 1);
    umag(i) = temp(end);
    umag_avg(i) = mean(u_mag(idx));
    
    [temp, I] = maxk(un_mag(idx), 1);
    unorm(i) = temp(end);
    unorm_avg(i) = mean(un_mag(idx));
    
    temp = un(idx);
    unorm_vec(i) = temp(I(end));
    unorm_vec_avg(i) = mean(temp);
    
    temp = maxk(ut_mag(idx), 1);
    utang(i) = temp(end);
    utang_avg(i) = mean(ut_mag(idx));
    
end

data.dist = dist;
data.umag = umag;
data.unorm = unorm;
data.utang = utang;
data.unorm_vect = unorm_vec;

data.umag_avg = umag_avg;
data.unorm_avg = unorm_avg;
data.utang_avg = utang_avg;
data.unorm_vect_avg = unorm_vec_avg;

[temp, I] = maxk(u_mag, 20);
data.max_umag = temp;
data.max_d = d(I);
data.max_unorm = un_mag(I);
data.max_utang = ut_mag(I);
end


function [pts, d, norm_vec] = generate_rand_pts(min_xyz, max_xyz, density, fv, L)

% Computer points all in some range of cell surface to reduce the
% computation time. So we will find points within 70 micron of the cell
% surface
cell_min = min(fv.vertices) - 30;
cell_max = max(fv.vertices) + 30;


% Find the true min_xyz taking into account the prvious image min_xyz and
% similar for max_xyz
min_xyz = max(min_xyz, cell_min);
max_xyz = min(max_xyz, cell_max);

% Number of points to generate
cube_edges_length = max_xyz - min_xyz;
nPts = round(density*prod((max(cube_edges_length))^3));
pts = rand(nPts, 3)*max(cube_edges_length);
idx = all(pts < cube_edges_length, 2);
pts = pts(idx, :) + min_xyz;

% Remove points from inside the cell
% tic
% in_pts_idx = inpolyhedron(fv, pts, 'FlipNormals',true);
% toc
% pts = pts(~in_pts_idx, :);

% % % Distance of each point to the cell surface
[d,surf_pts] = point2trimesh(fv, 'QueryPoints', pts);


% Save the points which are outside the cell surface and at most 70 um away
idx = d>0 & d<30;
pts = pts(idx, :);
d = d(idx);
surf_pts = surf_pts(idx, :);

% Compute normal vector to this point
norm_vec = pts - surf_pts;
norm_vec = norm_vec./(sqrt(sum(norm_vec.^2, 2))); %Unit vector

% Sort array by distance
[d, idx] = sort(d);
pts = pts(idx, :);
surf_pts = surf_pts(idx, :);
end


function [temp] = plot_coneplot(x, u, fv)

umag = sqrt(sum(u.^2,2));

% save cell surface
temp{1} = cell(4,1);
temp{1}{1,1} = fv.faces-1;
temp{1}{2,1} = fv.vertices;

% save displacement
temp{2} = cell(2,1);
temp{2}{1,1} = u;
temp{2}{2,1} = x;

%save max displacment
temp{3} = max(umag);

% save cell_center
temp{4} = mean(fv.vertices);


fv.vertices = fv.vertices(:,[2,1,3]); % convert to xyz format

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
% colorbar
% caxis([0 max(cmax)])
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
set(gca,'color','none','Visible','off')
set(h,'color','k');
close all;
end
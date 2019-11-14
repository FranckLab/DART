function [pts] = compute_disphist_dvc(fv, L, cell_idx, u, calib)
%Compute displacement histogram from the cell surface.


% Parameters
min_volume = 100;
density = 0.025;


% Load displacement
m = u.m;
u = u.u;
[u{[2,1,3]}] = u{[1,2,3]};
x0 = u.x0(u.track>0, :);
x1 = u.x1(u.track(u.track>0), :);
u = x1-x0;
x = num2cell(x1, 1);
% Remove rigid motion and calibrate displacement
u = removeRigidDriftRotation(u, x);

% Calibrate into proper units
fv.vertices = fv.vertices.*calib;
fv.faces = fv.faces(:,[1,3,2]);% flip faces
u = u.*calib;
x = num2cell(x1.*calib, 1);

% Compute cell volume and area
cell_volume = triVolume(fv.faces, fv.vertices);
cell_area = triArea(fv.faces, fv.vertices);



if cell_volume > min_volume
    
    % Generate random points around the cell in the valid range
    min_xyz = [30,30,15].*calib;
    max_xyz = (size(L) - [30,30,20]).*calib;
    [pts, d, n_vec] = generate_rand_pts(min_xyz, max_xyz, density, fv, L==cell_idx);
    
    for j = 1:3
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
    plot_coneplot(pts, u_pts, fv)
    
    % Compute histogram
    bins = 0:3:30
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
ut_mag = sqrt(sum(u_tang.^2 ,2));


% Compute histogram
umag = zeros(length(bins)-1, 1);
unorm = zeros(length(bins)-1, 1);
utang = zeros(length(bins)-1, 1);
dist = zeros(length(bins)-1, 1);

for i = 1:length(bins)-1
    dist(i) = mean(bins(i:i+1));
    idx = d>=bins(i) & d<bins(i+1);
    
    temp = maxk(u_mag(idx), 3);
    umag(i) = temp(end);
    
    temp = maxk(un_mag(idx), 3);
    unorm(i) = temp(end);
    
    temp = maxk(ut_mag(idx), 3);
    utang(i) = temp(end);

end

data.umag = umag;
data.unorm = unorm;
data.utang = utang;
data.dist = dist;
end


function [pts, d, norm_vec] = generate_rand_pts(min_xyz, max_xyz, density, fv, L)

% Computer points all in some range of cell surface to reduce the
% computation time. So we will find points within 70 micron of the cell
% surface
cell_min = min(fv.vertices) - 40;
cell_max = max(fv.vertices) + 40;


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
tic
[d,surf_pts] = point2trimesh(fv, 'QueryPoints', pts);
toc

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


function [hc] = plot_coneplot(x, u, fv)

fv.vertices = fv.vertices(:,[2,1,3]); % convert to xyz format
umag = sqrt(sum(u.^2,2));
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
end
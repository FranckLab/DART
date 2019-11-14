function [fn, c, azimuth, elevation] = calculateFaceNormals(tr)

fn = faceNormals(tr);
c = incenters(tr);
nx = fn(:,1);
ny = fn(:,2);
nz = fn(:,3);
% [azimuth,elevation,~] = cart2sph(nx,ny,nz);
nx = nx + eps;
azimuth = atan(ny./nx);
elevation = atan(nz./sqrt(nx.^2 + ny.^2));

azimuth = azimuth * 180/pi;
elevation = elevation * 180/pi;
end
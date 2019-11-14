
clear;
R = [10,20,30];
[x0, y0, z0] = ellipsoid(20,20,20,10,20,10,50);

% [f, v] = ellipse_tri('ico',4,[10,20,10]);

R = vrrotvec2mat([1 1 1, 45*pi/180]);

for i = 1:size(x0,1)
    for j = 1:size(x0,2)
        x_(i,j,:) = [x0(i,j), y0(i,j), z0(i,j)] * R;
        
    end
end
x = x_(:,:,1); y = x_(:,:,2); z = x_(:,:,3);
[v, f4] = surfToMesh(x,y,z);
f = triangulateFaces(f4);

figure;
trisurf(f,v(:,1),v(:,2),v(:,3)); axis image;

info = triInfo(f,v);

info.maxChordLengh*2






axis image;
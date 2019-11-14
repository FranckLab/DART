% u{1} = repmat(reshape(1:128,1,1,[]),[256 512 1]);
[m{1} m{2} m{3}] = meshgrid(1:256,1:256,1:256);
OS = [128,128,128];
u = sqrt(   (m{1} - OS(1)).^2 + (m{2} - OS(2)).^2 + (m{3} - OS(3)).^2   );

sizeI = size(u);

C = [64 64 64];
[f,v] = sphere_tri('ico',4,32);
v = bsxfun(@plus, v, C);

[n, ~, p] =  triNormals(f,v);
 
 uTri = triInterp(p,u);

 %%
figure;

hold on;
h1 = trisurf(f,v(:,1),v(:,2),v(:,3),uTri');
sliceIdx = sizeI/2; 
h2  = slice(m{1},m{2},m{3},u,  sliceIdx(1),sliceIdx(2),sliceIdx(3) );
alpha(h2, 0.5)       %# note: this switches to OpenGL renderer

axis equal
xlabel('x');
ylabel('y');
zlabel('z');
colorbar

set(h1,'linestyle','none');
set(h2,'linestyle','none');

 hold off
 view(3)
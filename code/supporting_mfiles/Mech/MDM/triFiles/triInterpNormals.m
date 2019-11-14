function [fn, s] = triInterpNormals(tri,s,plot_)

maxTime = length(tri);
fn = cell(maxTime,1); 
for t = 1:maxTime
    disp(['Interpolating normals onto displacement grid (', num2str(t) ,'/', num2str(maxTime),')']);
    [fn_, s_] = funInterpNormals(tri{t},s(t,:,:),plot_);
    fn{t,1} = fn_{1}; fn{t,2} = fn_{2}; fn{t,3} = fn_{3};
    s{t,3} = s_{3};
end

    
end

function [fn, s] = funInterpNormals(tri,s,plot_)
% tr:  delauney tiranguluzation of the surface. see TriRep()
% s{1}, s{2}: meshgrid that the tractions are defined on
% fnX, fnY, fnZ: interpolated surface
% inCenter: interpolated center positions of the triRep

triCenters = incenters(tri);
fn_ = faceNormals(tri);


% trInterp = TriScatteredInterp(mX,mY,mZ,'natural');
% trFunction = TriScatteredInterp(trSmooth.X(:,1),trSmooth.X(:,2),trSmooth.X(:,3),'natural');
centersFunction = TriScatteredInterp(triCenters(:,1),triCenters(:,2),triCenters(:,3));
fnXFunction = TriScatteredInterp(triCenters(:,1),triCenters(:,2),fn_(:,1));
fnYFunction = TriScatteredInterp(triCenters(:,1),triCenters(:,2),fn_(:,2));
fnZFunction = TriScatteredInterp(triCenters(:,1),triCenters(:,2),fn_(:,3));


% [mXInt,mYInt] = meshgrid(1:8:sizeI(1),1:8:sizeI(2));
% mZInt = trFunction(mXInt,mYInt);
s{3} = centersFunction(s{1},s{2});
fn{1} = fnXFunction(s{1},s{2});
fn{2} = fnYFunction(s{1},s{2});
fn{3} = fnZFunction(s{1},s{2});

s{3} = inpaint_nans(s{3});
fn{1} = inpaint_nans(fn{1});
fn{2} = inpaint_nans(fn{2});
fn{3} = inpaint_nans(fn{3});

if plot_ == 1
figure
trisurf(tri); shading flat; 
hold on;
quiver3(s{1},s{2},s{3},fn{1},fn{2},fn{3});
view(3);
axis equal
hold off;
colorbar
end


end



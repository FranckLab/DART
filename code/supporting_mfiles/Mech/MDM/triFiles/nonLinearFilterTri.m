 function tr_ = nonLinearFilterTri(tr, dim, iterations, fun, plot_)
tri = tr.ConnectivityList;
tr_ = tr;
nV = size(tr.Points,1);
vA = vertexAttachments(tr_,(1:nV)');

%% get Neighbors of vertices
nbr = cell(1,nV);
for i = 1:nV
    nbr{i} = (tr_(vA{i},:));
end
nbr = cellfun(@unique,nbr,'UniformOutput',0);

%% set function for each neigbors
stdPlot = zeros(1,iterations);
for i = 1:iterations
    
    msg = ['Applying median filter to mesh (Iteration: ', num2str(i),'/', num2str(iterations),')'];
    disp(msg);
    x = tr_.Points(:,1); y = tr_.Points(:,2); z = tr_.Points(:,3);
    
  
    
   
    if dim(1) == 1, x = setFun(x,nbr,fun); end
    if dim(2) == 1, y = setFun(y,nbr,fun); end
    if dim(3) == 1, z = setFun(z,nbr,fun); end
    
    tr_ = triangulation(tri,[x y z]);
    stdPlot(i) = mean2(std(tr_.Points,1));
end

if plot_ == 1
    figure;
subplot(2,2,1);
trisurf(tr_); shading flat
axis equal; colorbar;
title(msg);    
    
subplot(2,2,2)
plot(stdPlot/max(stdPlot(:)));
title('Convergence');
xlabel('Iteration');
ylabel('Normalized standard deviation of meshgrid coordinates');
end

 end

function coord = setFun(coord,nbr,fun)
coord_ = subCellRef(coord,nbr);
coord = cellfun(fun,coord_);
end

function p = subCellRef(v,ref)

nP = length(v);
p = cell(nP,1);
for i = 1:nP
    p{i} = v(ref{i});
end
end
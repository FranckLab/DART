function ITri = triInterp(V,I,m,method)
if nargin < 4, method = 'cubic'; end
if nargin < 3, m = 1; end



if ~iscell(I), I = {I}; end

sizeI = size(I{1});

if iscell(m)
    if ismatrix(m{1})
       mIdx = m; 
       [m{1}, m{2}, m{3}] = ndgrid(mIdx{:});
    end
else
    dm = m; clear m;
    if numel(dm) == 1, dm = dm*[1 1 1]; end
    mIdx = cell(1,3);
    for i = [1 2 3], mIdx{i} = dm(i):dm(i):dm(i)*sizeI(i); end    
    [m{1}, m{2}, m{3}] = ndgrid(mIdx{:});
end



ITri = zeros(size(I,1),size(I,2), size(V,1));

for i = 1:size(I,1)
    for j = 1:size(I,2)
        interpolant = griddedInterpolant(m{1}, m{2}, m{3}, I{i,j}, method);
        ITri(i,j,:) = interpolant(V(:,2),V(:,1),V(:,3));
    end
end

ITri = squeeze(ITri);
if ismatrix(ITri), ITri = ITri'; end


end
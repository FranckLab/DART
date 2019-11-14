function u = changeUnits(u,conversion)
if ~iscell(u), u = {u}; end

if length(conversion) == 1, conversion = conversion*ones(1,numel(u{1})); end


for i = 1:length(u)
    nDim = length(conversion);
    u{i} = u{i}(1:nDim);
    for j = 1:nDim
         u{i}{j} = u{i}{j}*conversion(j);  
    end
end


end
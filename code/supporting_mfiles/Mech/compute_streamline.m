function [xyz, u_xyz] = compute_streamline(ugrid, XYZ, spts, BW, xy, w)
%Compute streamline

% idx = randperm(length(spts),length(spts));
% spts = spts(idx, :); clear idx

xyz = cell(size(spts,1),1);
for i = 1:size(spts, 1)
    xyz{i} = funstream3(ugrid, XYZ, spts(i,:));
end
xyz = cell2mat(xyz);
idx = 1:24:length(xyz);
xyz = xyz(idx, :);

% if size(xyz, 1) > 10000
%     xyz = xyz(1:10000,:);
% end

% xyz_round = round(xyz) + 2*[w,w,0] - xy;


for j = 1:3
    u_xyz(:,j) = interp3(XYZ{:}, ugrid{j}, xyz(:,1), xyz(:,2), xyz(:,3));
end

end


function xyz = funstream3(u, XYZ, spts)
%Compute streamline
% [u{:}] = u{[2,1,3]};
% Combine streamline in both the directions
xyz1 = stream3(XYZ{:}, u{:}, spts(:,1), spts(:,2), spts(:,3), 0.2);
xyz2 = stream3(XYZ{:}, -u{1}, -u{2}, -u{3}, spts(:,1), spts(:,2), ...
                                                        spts(:,3), 0.2);
xyz = [xyz1{1}; xyz2{1}];

end

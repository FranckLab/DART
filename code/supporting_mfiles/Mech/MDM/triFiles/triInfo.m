function info = triInfo(f,v)

info = struct;

v0 = v(f(:,1),:);
v1 = v(f(:,2),:);
v2 = v(f(:,3),:);

Ni = cross(v1 - v0,v2 - v0);  % normal vector (not normalized)
Ai = 0.5*sqrt(sum(Ni.^2,2)); % area for each triangle
A = sum(Ai); % total surface area
Ci = bsxfun(@times, (v0+v1+v2)/3,Ai); % centroid for each triangle
C = sum(Ci)/A;  % centroid

Vi = 1/6*abs(sum(Ni.*v0,2)); % volume for each pyramid
V = sum(Vi);  % total volume

% calculate max chord length
nChords = size(v,1);
v_ = bsxfun(@minus, v, mean(v));
v_(abs(v_) < eps(1)) = 0;

maxChordLength_ = zeros(nChords,1);

for i = 1:nChords
    chords = -bsxfun(@minus,v_,v_(i,:));
    chordMag = sqrt(sum(chords.^2,2));
    maxChordLength_(i) = max(chordMag);
end

maxChordLength = max(chordMag);

% calculate major,minor axis
% FOR ANOTHER TIME
% L = sqrt(sum(v(diametricIdx,:).^2,2));
% majorAxis = max(L);
% minorAxis = min(L);


info.area = A;
info.volume = V;
info.centroid = C;
% info.distance = L;
% info.minorAxis = minorAxis;
% info.majorAxis = majorAxis;
info.normals = bsxfun(@rdivide,Ni,sqrt(sum(Ni.^2,2)));
info.maxChordLength = maxChordLength;

% calculate mean radius
info.meanRadius = mean(mean(abs(v_)));

end



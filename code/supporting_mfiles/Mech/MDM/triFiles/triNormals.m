function [FN, VN, P] = triNormals(f,v)

if nargin == 1, 
    if isstruct(f) % if input is a structure from isosurface
     F = f.faces;
     V = f.vertices;
     TR = triangulation(F,V);
    else % if input is a triangulation
        TR = f;
    end
else
    TR = triangulation(f,v);
end

FN = faceNormal(TR);
VN = vertexNormal(TR);
P = incenter(TR);

% flip normal directions to inward for matrix calculations
% FN = -FN;
% VN = -VN;
end

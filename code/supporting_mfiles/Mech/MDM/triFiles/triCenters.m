function P = triCenters(f,v)

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

P = incenter(TR);
end

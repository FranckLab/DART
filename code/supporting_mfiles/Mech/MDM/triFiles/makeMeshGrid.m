function [msg,nx,ny,nz] = makeMeshGrid(x,y,z,data)
%XYZVCHECK  Check arguments to 3D scalar data routines.
%   [MSG,X,Y,Z] = XYZVCHECK(X,Y,Z,V) checks the input arguments
%   and returns either an error message structure in MSG or valid
%   X,Y,Z. The ERROR function describes the format and use of the
%   error structure.
%
%   See also ERROR

%   Copyright 1984-2011 The MathWorks, Inc. 
%   $Revision: 1.6.4.3 $  $Date: 2011/08/13 17:30:46 $

msg = struct([]);
nx = x;
ny = y;
nz = z;

sz = size(data);

if ndims(data)~=3
  msg(1).identifier = 'MATLAB:xyzvcheck:VNot3D';
   msg(1).message = getString(message(msg(1).identifier));
  return
end
if min(sz)<2
  msg(1).identifier = 'MATLAB:xyzvcheck:VPlanar';
   msg(1).message = getString(message(msg(1).identifier)); 
  return
end

nonempty = ~[isempty(x) isempty(y) isempty(z)];
if any(nonempty) && ~all(nonempty)
  msg(1).identifier = 'MATLAB:xyzvcheck:XYZMixedEmpty';
   msg(1).message = getString(message(msg(1).identifier));
  return;
end

if ~isempty(nx) && ~isequal(size(nx), sz)
  nx = nx(:);
  if length(nx)~=sz(2)
    msg(1).identifier = 'MATLAB:xyzvcheck:XVSizeMismatch';
   msg(1).message = getString(message(msg(1).identifier));
    return
  else
    nx = repmat(nx',[sz(1) 1 sz(3)]);
  end
end

if ~isempty(ny) && ~isequal(size(ny), sz)
  ny = ny(:);
  if length(ny)~=sz(1)
    msg(1).identifier = 'MATLAB:xyzvcheck:YVSizeMismatch';
   msg(1).message = getString(message(msg(1).identifier));
    return
  else
    ny = repmat(ny,[1 sz(2) sz(3)]);
  end
end

if ~isempty(nz) && ~isequal(size(nz), sz)
  nz = nz(:);
  if length(nz)~=sz(3)
    msg(1).identifier = 'MATLAB:xyzvcheck:ZVsizeMismatch';
   msg(1).message = getString(message(msg(1).identifier));
    return
  else
    nz = repmat(reshape(nz,[1 1 length(nz)]),[sz(1) sz(2) 1]);
  end
end


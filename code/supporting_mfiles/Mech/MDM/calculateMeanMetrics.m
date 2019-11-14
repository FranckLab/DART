function mechanics = calculateMeanMetrics(f,v,uTri)

if ~iscell(f), f = {f}; v = {v}; uTri = {uTri}; end

mechanics = struct;
for i = 1:length(uTri)
    mechanics = funCalculateMeanMetrics(f{i},v{i},uTri{i}, mechanics, i);
end


end

function mechanics = funCalculateMeanMetrics(f,v,uTri, mechanics, i)

[fn, fv, fp]  = triNormals(f,v);  % normal vector on faces

if all(size(uTri) == size(f)), % u defined on faces
    n = fn;
    p = fp; 
else % u defined on vertices
    n = fv;
    p = v;
end

% mean velocity of center of mass calculation

x = bsxfun(@minus, p, triCentroid(f,v));

xi_uj = vecDyadic(x,uTri);
xi_uj_nj = zeros(size(xi_uj,3),3);
for j = 1:size(xi_uj,3)
    xi_uj_nj(j,:) = xi_uj(:,:,j)*n(j,:)';
end

triV = triVolume(f,v);
triA = triArea(f,v);
uCM = 1/triV*triIntegrate(f,v, xi_uj_nj);

ui_nj = vecDyadic(uTri,n);
I = eye(3,3);
F = I + (1/triV)*triIntegrate(f,v, ui_nj);
J = det(F);
C = F'*F;
B = F*F';

if all(~isnan(C(:)))
    U = sqrtm(C);
    V = sqrtm(B);
else
    U = nan(3,3);
    V = nan(3,3);
end

E = 1/2*(C - I);
e = triIntegrate(f,v, ui_nj);
e = 1/(2*triV)*(e + e');
R = V\F;
axisAngleR = vrrotmat2vec(R);

if all(~isnan(E(:)))
    [eigVecE, eigValE] = eig(E);
    eigValE = diag(eigValE);
else
    eigValE = nan(1,3);
    eigVecE = nan(3,3);
end

mechanics.triV(i,1) = triV;
mechanics.triA(i,1) = triA;
mechanics.F(:,:,i) = F;
mechanics.J(i,1) = J;
mechanics.U(:,:,i) = U;
mechanics.C(:,:,i) = C;
mechanics.E(:,:,i) = E;
mechanics.eigValE(i,:) = eigValE;
mechanics.eigValENorm(i,:) = eigValE/norm(eigValE);
mechanics.eigVecE(:,:,i) = eigVecE;
mechanics.e(:,:,i) = e;
mechanics.R(:,:,i) = R;
mechanics.axisAngleR(i,:) = axisAngleR;
mechanics.uCM(i,:) = uCM;



end
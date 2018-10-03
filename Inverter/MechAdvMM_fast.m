function [MA, dMA] = MechAdvMM_fast(x, Evec, F, constraints, ind_out, Ex, Ey)

% Function evaluates the mechanical advantage

global  p_SIMP u_def W_filt xn xm Em beta

nu = 0.3;nelx=Ex;nely=Ey;nelem=Ex*Ey;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);


Hs=ones(nelem,1);H=W_filt;
x=reshape(x,nelem,xm);
xFilt=(H*x)./(Hs*ones(1,xm));
[xHsvi, dxHsvi] = Heaviside(xFilt(:), beta);
[rMM, drMM] = r_func_fast(xHsvi, nelem, xm, Em, p_SIMP);
drMM = drMM.*repmat(reshape(dxHsvi,nelem,xm),1,Em+1);
Emm = rMM*Evec';
dEmm=zeros(nelem,xm);
for i=1:xm
    dEmm(:,i) = drMM(:,i:xm:end) * Evec';
end

alldofs = [1:2*(nely+1)*(nelx+1)];
fixeddofs = union(constraints, ind_out);
freedofs = setdiff(alldofs,fixeddofs);
U = zeros(2*(nely+1)*(nelx+1),1);


sK = reshape(KE(:)*(Emm'),64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
U(freedofs) = K(freedofs,freedofs)\F(freedofs,:);
U(fixeddofs, 1) = 0;
F_int = K*U;
F_out = -F_int(ind_out);
F_mag = sum(F);
MA = F_out/F_mag;


psi = zeros(2*(nely+1)*(nelx+1), 1);
psi(ind_out) = -1/F_mag;
psi(freedofs) = -K(freedofs, freedofs)\(K(freedofs, fixeddofs)*psi(fixeddofs));



MA_ce = reshape(sum((psi(edofMat)*KE).*U(edofMat),2),nelem,1);
dMA = dEmm.*repmat(MA_ce,1,xm);
dMA=([dMA./(Hs*ones(1,xm))]'*H)';
dMA=reshape(dMA,nelem*xm,1);





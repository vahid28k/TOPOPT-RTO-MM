% This Matlab code was written by Vahid Keshavarzzadeh and Kai A. James
function [Q, dQ]=comp_fast(x,F,fixeddofs,lambda)
%% MATERIAL PROPERTIES
global  p_SIMP Em beta Ex Ey xm Evec W_filt;

nelx=Ex;
nely=Ey;

E0 = 1;
Emin = 1e-9;
nu = 0.3;
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
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
%F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
%fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);

%% Density Process
%x = repmat(volfrac,nely,nelx);


nelem=Ex*Ey;
Hs=ones(nelem,1);H=W_filt;
x=reshape(x,nelem,xm);
xFilt=(H*x)./(Hs*ones(1,xm));
[xHsvi, dxHsvi] = Heaviside(xFilt(:), beta);
[rMM, drMM] = r_func_fast(xHsvi, nelem, xm, Em, p_SIMP);
drMM = drMM.*repmat(reshape(dxHsvi,nelem,xm),1,Em+1);

load XW6.mat XW;
comp=zeros(1,size(XW,1));dcomp=zeros(nelem*xm,size(XW,1));
for i_q=1:size(XW,1)
    
Evec_quad=exp(log(Evec)+(0.05*log(Evec).*[0 XW(i_q,1:Em)]));
    
Emm = rMM*Evec_quad';
dEmm=zeros(nelem,xm);
for i=1:xm
    dEmm(:,i) = drMM(:,i:xm:end) * Evec_quad';
end
%% FE-ANALYSIS
sK = reshape(KE(:)*(Emm'),64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
U(freedofs) = K(freedofs,freedofs)\F(freedofs);
%% OBJECTIVE FUNCTION AND SENSITIVITY 
ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nelem,1);
c = sum(Emm.*ce);
dc = -dEmm.*repmat(ce,1,xm);
dc=([dc./(Hs*ones(1,xm))]'*H)';
dc=reshape(dc,nelem*xm,1);

comp(1,i_q)=c;
dcomp(:,i_q)=dc;

end

cmean=comp*XW(:,4);cmean2=comp.^2*XW(:,4);
dcmean=dcomp*XW(:,4);cdcmean=(repmat(comp,nelem*xm,1).*dcomp)*XW(:,4);

%lambda weight of standard deviation in Robust Design i.e. obj function: E[C]+ \lambda \sigma(C)
Q=cmean+lambda*sqrt(cmean2-cmean^2);dQ=dcmean+lambda*([cdcmean-cmean*dcmean]/sqrt(cmean2-cmean^2));

%sigm=sqrt(cmean2-cmean^2);

 




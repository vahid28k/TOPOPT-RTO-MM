%%%% This Matlab code was written by Vahid Keshavarzzadeh and Kai A. James
function [Q1, dQ1, Q2,dQ2, Q, dQ]=comp_fast(x,F1,F2,fixeddofs,ind_out,lambda)
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
comp1=zeros(1,size(XW,1));dcomp1=zeros(nelem*xm,size(XW,1));
comp2=zeros(1,size(XW,1));dcomp2=zeros(nelem*xm,size(XW,1));
madv=zeros(1,size(XW,1));dmadv=zeros(nelem*xm,size(XW,1));
for i_q=1:size(XW,1)
    
Evec_quad=exp(log(Evec)+(0.02*log(Evec).*[0 XW(i_q,1:Em)]));

Emm = rMM*Evec_quad';
dEmm=zeros(nelem,xm);
for i=1:xm
    dEmm(:,i) = drMM(:,i:xm:end) * Evec_quad';
end
%% FE-ANALYSIS
sK = reshape(KE(:)*(Emm'),64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2;

U1 = zeros(2*(nely+1)*(nelx+1),1);
U2 = zeros(2*(nely+1)*(nelx+1),1);
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
U1(freedofs) = K(freedofs,freedofs)\F1(freedofs);
U2(freedofs) = K(freedofs,freedofs)\F2(freedofs);
%% OBJECTIVE FUNCTION AND SENSITIVITY 
ce1 = reshape(sum((U1(edofMat)*KE).*U1(edofMat),2),nelem,1);
c1 = sum(Emm.*ce1);
dc1 = -dEmm.*repmat(ce1,1,xm);
dc1=([dc1./(Hs*ones(1,xm))]'*H)';
dc1=reshape(dc1,nelem*xm,1);
ce2 = reshape(sum((U2(edofMat)*KE).*U2(edofMat),2),nelem,1);
c2 = sum(Emm.*ce2);
dc2 = -dEmm.*repmat(ce2,1,xm);
dc2=([dc2./(Hs*ones(1,xm))]'*H)';
dc2=reshape(dc2,nelem*xm,1);

U = zeros(2*(nely+1)*(nelx+1),1);
alldofs = [1:2*(nely+1)*(nelx+1)];
fixeddofsMA = union(fixeddofs, ind_out);
freedofs = setdiff(alldofs,fixeddofsMA);
U(freedofs) = K(freedofs,freedofs)\F1(freedofs,:);
U(fixeddofsMA, 1) = 0;
F_int = K*U;
F_out = -F_int(ind_out);
F_mag = sum(F1);
MA = F_out/F_mag;
psi = zeros(2*(nely+1)*(nelx+1), 1);
psi(ind_out) = -1/F_mag;
psi(freedofs) = -K(freedofs, freedofs)\(K(freedofs, fixeddofsMA)*psi(fixeddofsMA));
MA_ce = reshape(sum((psi(edofMat)*KE).*U(edofMat),2),nelem,1);
dMA = dEmm.*repmat(MA_ce,1,xm);
dMA=([dMA./(Hs*ones(1,xm))]'*H)';
dMA=reshape(dMA,nelem*xm,1);

comp1(1,i_q)=c1;
dcomp1(:,i_q)=dc1;
comp2(1,i_q)=c2;
dcomp2(:,i_q)=dc2;
madv(1,i_q)=MA;
dmadv(:,i_q)=dMA;

end

cmean=madv*XW(:,4);cmean2=madv.^2*XW(:,4);
dcmean=dmadv*XW(:,4);cdcmean=(repmat(madv,nelem*xm,1).*dmadv)*XW(:,4);

%lambda weight of standard deviation in Robust Design i.e. obj function: E[C]+ \lambda \sigma(C)
Q=cmean+0*sqrt(cmean2-cmean^2);dQ=dcmean+0*([cdcmean-cmean*dcmean]/sqrt(cmean2-cmean^2));

c1mean=comp1*XW(:,4);c1mean2=comp1.^2*XW(:,4);
dc1mean=dcomp1*XW(:,4);cdc1mean=(repmat(comp1,nelem*xm,1).*dcomp1)*XW(:,4);
Q1=c1mean+lambda*sqrt(c1mean2-c1mean^2);dQ1=dc1mean+lambda*([cdc1mean-c1mean*dc1mean]/sqrt(c1mean2-c1mean^2));


c2mean=comp2*XW(:,4);c2mean2=comp2.^2*XW(:,4);
dc2mean=dcomp2*XW(:,4);cdc2mean=(repmat(comp2,nelem*xm,1).*dcomp2)*XW(:,4);
Q2=c2mean+lambda*sqrt(c2mean2-c2mean^2);dQ2=dc2mean+lambda*([cdc2mean-c2mean*dc2mean]/sqrt(c2mean2-c2mean^2));
 
 





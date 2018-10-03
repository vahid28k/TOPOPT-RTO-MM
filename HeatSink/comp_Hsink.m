% This Matlab code was written by Vahid Keshavarzzadeh and Kai A. James
function [Q, dQ]=comp_Hsink(x,F,fixeddofs, edofMat,lambda)
%% MATERIAL PROPERTIES
global  p_SIMP Em beta Ex Ey xm Evec W_filt;

nelx=Ex;
nely=Ey;

E0 = 1;
Emin = 1e-3;

KE =[2/3 -1/6 -1/3 -1/6;
    -1/6 2/3 -1/6 -1/3;
    -1/3 -1/6 2/3 -1/6;
    -1/6 -1/3 -1/6 2/3]; 
iK = reshape(kron(edofMat,ones(4,1))',16*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,4))',16*nelx*nely,1);

U = zeros((nely+1)*(nelx+1),1);
alldofs     = [1:(nely+1)*(nelx+1)];
freedofs    = setdiff(alldofs,fixeddofs);


%% Density Process


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
    sK = reshape(KE(:)*(Emm'),16*nelx*nely,1);
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
%csigma=sqrt(cmean2-cmean^2);






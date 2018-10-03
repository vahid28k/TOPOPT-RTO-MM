function [Q, dQ] = volRectMM_fast(x, xm, Ex, Ey, lx, ly)

global W_filt Evec xn Em beta

% function computes volume fraction of a rectangualr domain with uniform
% mesh

ne = Ex*Ey;
V = 0;
dV = zeros(ne*xm, 1);


x=reshape(x,ne,xm);
xFilt=(W_filt*x);
[xHsvi, dxHsvi] = Heaviside(xFilt(:), beta);
[rMM, drMM] = r_func_fast(xHsvi, ne, xm, Em, 1);
drMM = drMM.*repmat(reshape(dxHsvi,ne,xm),1,Em+1);




%%%%%%%%%%%% COST FUNCTION %%%%%%%%%%%%%%

Emax = max(Evec(2:Em+1));
%[rMM, drx] = r_func(xf, ne, xm, Em, 1);
%rMM = W_filt*rMM;

load XW6.mat XW;
vol=zeros(1,size(XW,1));dvol=zeros(ne*xm,size(XW,1));
for i_q=1:size(XW,1)
    Evec_quad=exp(log(Evec)+(0.05*log(Evec).*[0 XW(i_q,1:Em)]));
    
    V = (rMM*Evec_quad')/(Emax*ne);V=sum(V);
    dVmm=zeros(ne,xm);
    for i=1:xm
        dVmm(:,i) = (drMM(:,i:xm:end) * Evec_quad')/(Emax*ne);
    end
    dVmm=(dVmm'*W_filt)';
    dV=reshape(dVmm,ne*xm,1);
    
    vol(1,i_q)=V;
    dvol(:,i_q)=dV;
end
Q=vol*XW(:,4);
dQ=dvol*XW(:,4);



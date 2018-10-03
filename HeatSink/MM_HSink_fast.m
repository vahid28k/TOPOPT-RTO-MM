% Main script for optimizing a mechanical inverter with multiple materials
% Heaviside version
clc;clear all;close all;

global  p_SIMP W_filt constraints nodeMap Ex Ey lx ly density rho_min beta xval u_def xn xm Em Evec

% Set material parameters
% 
% Set domain parameters

Ex = 80; % number of elements in X direction
Ey = 80; % number of elements in Y direction

lx = 80.0; % domain length in X direction
ly = 80.0; % domain length in Y direction

nel = 2; % linear elements (used in defining BCs)
ne = Ex*Ey; % number of elements
beta = 0.0; % Heaviside regularization constant

p_SIMP = 3.0; % SIMP penalization parameter
rho_min = 1e-3;
density = 1.0;% volumetric mass density of material

xn = ne; % number of elements
xm = 2;  % number of design variables per element
Em = 3;  % number of used/design materials (not including void)

Evec = [1e0, 1e8, 2e8, 3e8]; % Young's moduli of candidate materials (first number is modulus for void regions)
v = 0.3; % Poisson's ratio

x0 = 0.4*ones(xn*xm, 1); % initial uniform density value (volume fraction)
%x0 = random('unif', 0.0, 1.0, xn*xm, 1); % random initial distribution for testing sensitivities


% Initialize boundary conditions
   
F_ext = sparse((Ey+1)*(Ex+1),1);
F_ext(:,1) = 0.01;
bnc=4;
fixeddofs   = [Ey/2+1-bnc:2:Ey/2+1+bnc];
constraints=fixeddofs;


edofMat=[];
for elx = 1:Ex
    for ely = 1:Ey
        n1 = (Ey+1)*(elx-1)+ely;
        n2 = (Ey+1)* elx   +ely;
        edofMat =[edofMat; [n1 n2 n2+1 n1+1]];
    end
end


% Initialize density filter
nh = 3; % filter neighborhood (approximate radius measured in elements across)
W_filt = density_filter_fast(Ex, Ey, nh); % linear filter

xval = x0; % initialize design variables



%rho = W_filt*xval;
%[rhoH, drh] = Heaviside(rho, beta);

% Solve linear FEM problem; evaluate objective and sensitivities
lambda=1;
[C, dC] = comp_Hsink(xval, F_ext, constraints, edofMat,lambda);
[V, dV] = volRectMM_fast(xval, xm, Ex, Ey, lx, ly);
[B, dB] = binaryCon(xval, xm*xn);

obj_hist(1) = C;

C0 = C;
C_scale = 8/C0;
V_bound = 0.30; % volume constraint (global volume fraction)



% Initialize parameters used in MMA problem

m = 1; % 3 constraints
n = Ex*Ey;
epsimin = 1e-7;
xold1   = xval;
xold2   = xval;
xmin    = 0*ones(xn*xm, 1);
xmax    = ones(xn*xm, 1);
low     = xmin;
upp     = xmax;
c_mma   = 1000*ones(m, 1);
d_mma   = 1*ones(m, 1);
a0_mma  = 1;
a_mma   = 0*ones(m, 1);
outeriter = 0;
maxoutit  = 150;

kkttol  = 5.0e-4;
kktnorm = kkttol + 1;

f0val = C*C_scale;
df0dx = dC*C_scale;
fval = V - V_bound;
dfdx = dV;

%disp([' ************************************* ']);
disp([' outiter: ' sprintf('%4i', outeriter) ' kktnorm: ' sprintf('%6.4e', kktnorm) ' Compliance: ' sprintf('%10.7f', full(C)*C_scale) ...
      ' VolFrac: ' sprintf('%5.3f', V) ' Binary Coeff.: ' sprintf('%5.3f', B)]);
%disp([' ************************************* ']);

% Plot material distribution
xf=reshape(xval,ne,xm);
%xHS = Heaviside(xf, beta);
rMM = r_func_fast(xf, xn, xm, Em, p_SIMP);
Emm = rMM*Evec';
Emm = W_filt*Emm;


rhoBlock = reshape(Emm, Ey, Ex);
figure(1)
imagesc(real(rhoBlock));colormap(jet);axis equal;axis off;colorbar;
%snapshots(outeriter+1) = getframe;

% %*************
% % Hot start
% outeriter = 0;
% kktnorm = 1;
% %*************

while outeriter < maxoutit && kktnorm > kkttol

    outeriter = outeriter+1;
    
    %%%% The MMA subproblem is solved at the point xval:
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = mmasub(m,xn*xm,outeriter,xval,xmin,xmax,xold1,xold2,f0val,df0dx,fval',dfdx',low,upp,a0_mma,a_mma,c_mma,d_mma);
    
    %%%% Some vectors are updated:
    xold2 = xold1;
    xold1 = xval;
    xval  = xmma;
    
    
    % Solve linear FEM problem; evaluate objective and sensitivities
    
    [C, dC] = comp_Hsink(xval, F_ext, constraints, edofMat,lambda);
    [V, dV] = volRectMM_fast(xval, xm, Ex, Ey, lx, ly);
    [B, dB] = binaryCon(xval, xm*xn);
    
    obj_hist(outeriter+1) = C;
    
    f0val = C*C_scale;
    df0dx = dC*C_scale;
    fval = V - V_bound;
    dfdx = dV;

    %%%% The residual vector of the KKT conditions is calculated:
    [residu,kktnorm,residumax] = kktcheck(m,xn*xm,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,xmin,xmax,df0dx,fval',dfdx',a0_mma,a_mma,c_mma,d_mma);

    %disp([' ************************************* ']);
    disp([' outiter: ' sprintf('%4i', outeriter) ' kktnorm: ' sprintf('%6.4e', kktnorm) ' Compliance: ' sprintf('%10.7f', full(C)*C_scale) ...
          ' VolFrac: ' sprintf('%5.3f', V) ' Binary Coeff.: ' sprintf('%5.3f', B)]);
    %disp([' ************************************* ']);

    % Plot material distribution
    xf=reshape(xval,ne,xm);
    %xHS = Heaviside(xf, beta);
    rMM = r_func_fast(xf, xn, xm, Em, p_SIMP);
    Emm = rMM*Evec';
    Emm = W_filt*Emm;



    rhoBlock = reshape(Emm, Ey, Ex);
    figure(1)
    imagesc(real(rhoBlock));colormap(jet);axis equal;axis off;colorbar;
    %snapshots(outeriter+1) = getframe;
    
    % plot single material distribution
%     for j = 2:Em+1
%         rhoM = rMM(:, j);
%         rhoBlock = vec2array(rhoM, Ex, Ey);
%         rhoBlockSym(:, Ex+1:2*Ex) = rhoBlock;
%         for i = 1:Ex
%             rhoBlockSym(:, i) = rhoBlock(:, Ex+1-i);
%         end
% 
%         figure(j)
%         imagesc(real(-rhoBlockSym));colormap(gray);axis equal;axis off;colorbar;
%     end 

data_hist1(outeriter,:)=[full(C)*C_scale V];
data_hist2(:,outeriter)=xval;

end






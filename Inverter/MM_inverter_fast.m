% Main script for optimizing a mechanical inverter with multiple materials
 clc;clear all;close all;

global  p_SIMP W_filt constraints nodeMap Ex Ey lx ly density rho_min beta xval u_def xn xm Em Evec

% Set material parameters

% Set domain parameters

Ex = 48; % number of elements in X direction
Ey = 96; % number of elements in Y direction

lx = 2.5; % domain length in X direction
ly = 5.0; % domain length in Y direction

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

x0 = 0.5*ones(xn*xm, 1); % initial uniform density value (volume fraction)
%x0 = random('unif', 0.4, 0.8, xn*xm, 1); % random initial distribution for testing sensitivities

F_out = zeros(2*(Ex+1)*(Ey+1), 1);
F_ext = zeros(2*(Ex+1)*(Ey+1), 1);
ind_out=(Ey+1)*2;
F_ext(2,1)=0.5;F_ext(ind_out+2,1)=0.5;
F_out(ind_out,1) = 1;
constraints = union([1:2:2*(Ey+1)],[2*(Ex)*(Ey+1)+2]);



% Initialize density filter
nh =2.5; % filter neighborhood (approximate radius measured in elements across)
W_filt = density_filter_fast(Ex, Ey, nh); % linear filter


xval = x0; % initialize design variables

% Solve linear FEM problem; evaluate objective and sensitivities
lambda=0.1;
[C, dC, C2, dC2, MA, dMA] = comp_fast(xval, F_ext, F_out, constraints, ind_out,lambda);
[V, dV] = volRectMM_fast(xval, xm, Ex, Ey, lx, ly);
[B, dB] = binaryCon(xval, xm*xn);

MA0 = MA;
MA_scale = 1/MA0;

C0 = C;
C_scale = 1e6;

C20 = C2;
C2_scale = 1/C20;

C_bound = 0.75; % maximum allowable compliance (expressed as fraction of initial value)
C2_bound = 1; % maximum allowable compliance (expressed as fraction of initial value)
V_bound = 0.2; % volume constraint (global volume fraction)
B_bound = (1e-3); % uppoer bound on binary constraint

%*************************************************************************


% Initialize parameters used in MMA problem



%*************************************************************************


% Initialize parameters used in MMA problem

m = 3; % 3 constraints
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
maxoutit  = 200;

kkttol  = 5.0e-4;
kktnorm = kkttol + 1;

f0val = MA*MA_scale;
df0dx = dMA*MA_scale;
fval(1) = C*C_scale - C_bound;
dfdx = dC*C_scale;
fval(2) = C2*C2_scale - C2_bound;
dfdx(:, 2) = dC2*C2_scale;
fval(3) = V - V_bound;
dfdx(:, 3) = dV;
%fval(3) = B - B_bound;
%dfdx(:, 3) = dB;

%disp([' ************************************* ']);
disp([' outiter: ' sprintf('%4i', outeriter) ' kktnorm: ' sprintf('%6.4e', kktnorm) ' Mech adv: ' sprintf('%10.7f', MA) ' Comp constraint: ' sprintf('%10.7f', full(C)*C_scale/C_bound) ...
      ' Comp constraint 2: ' sprintf('%10.7f', full(C2)*C2_scale/C2_bound) ' VolFrac: ' sprintf('%5.3f', V) ' Binary Coeff.: ' sprintf('%5.3f', B)]);
%disp([' ************************************* ']);

% Plot material distribution

xf=reshape(xval,ne,xm);
rMM = r_func_fast(xf, xn, xm, Em, p_SIMP);
Emm = rMM*Evec';
Emm = W_filt*Emm;

rhoBlock = reshape(Emm, Ey, Ex);
    rhoBlockSym(:, Ex+1:2*Ex) = rhoBlock;
    for i = 1:Ex
        rhoBlockSym(:, i) = rhoBlock(:, Ex+1-i);
    end
% figure(1)
% imagesc(real(rhoBlockSym));colormap(jet);axis equal;axis off;colorbar;
%snapshots(outeriter+1) = getframe;

while outeriter < maxoutit && kktnorm > kkttol

    outeriter = outeriter+1;
    
    %%%% The MMA subproblem is solved at the point xval:
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = mmasub(m,xn*xm,outeriter,xval,xmin,xmax,xold1,xold2,f0val,df0dx,fval',dfdx',low,upp,a0_mma,a_mma,c_mma,d_mma);
    
    %%%% Some vectors are updated:
    xold2 = xold1;
    xold1 = xval;
    xval  = xmma;
    
    % Solve linear FEM problem; evaluate objective and sensitivities
    [C, dC, C2, dC2, MA, dMA] = comp_fast(xval, F_ext, F_out, constraints, ind_out,lambda);
    [V, dV] = volRectMM_vahid(xval, xm, Ex, Ey, lx, ly);
    [B, dB] = binaryCon(xval, xm*xn);
    
    f0val = MA*MA_scale;
    df0dx = dMA*MA_scale;
    fval(1) = C*C_scale - C_bound;
    dfdx = dC*C_scale;
    fval(2) = C2*C2_scale - C2_bound;
    dfdx(:, 2) = dC2*C2_scale;
    fval(3) = V - V_bound;
    dfdx(:, 3) = dV;
    %fval(3) = B - B_bound;
    %dfdx(:, 3) = dB;
    
    %%%% The residual vector of the KKT conditions is calculated:
    [residu,kktnorm,residumax] = kktcheck(m,xn*xm,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,xmin,xmax,df0dx,fval',dfdx',a0_mma,a_mma,c_mma,d_mma);

    %disp([' ************************************* ']);
    disp([' outiter: ' sprintf('%4i', outeriter) ' kktnorm: ' sprintf('%6.4e', kktnorm) ' Mech adv: ' sprintf('%10.7f', MA) ' Comp constraint: ' sprintf('%10.7f', full(C)*C_scale/C_bound) ...
          ' Comp constraint 2: ' sprintf('%10.7f', full(C2)*C2_scale/C2_bound) ' VolFrac: ' sprintf('%5.3f', V) ' Binary Coeff.: ' sprintf('%5.3f', B)]);
    %disp([' ************************************* ']);
    %ymma
    %zmma
    

    % Plot material distribution
    xf=reshape(xval,ne,xm);
    %xHS = Heaviside(xf, beta);
    rMM = r_func_fast(xf, xn, xm, Em, p_SIMP);
    Emm = rMM*Evec';
    Emm = W_filt*Emm;


    rhoBlock = reshape(Emm, Ey, Ex);
    rhoBlockSym(:, Ex+1:2*Ex) = rhoBlock;
    for i = 1:Ex
        rhoBlockSym(:, i) = rhoBlock(:, Ex+1-i);
    end
%     figure(1)
%     imagesc(real(rhoBlockSym));colormap(jet);axis equal;axis off;colorbar horizontal;
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

data_hist1(outeriter,:)=[full(C)*C_scale/C_bound full(C2)*C2_scale/C2_bound V];
data_hist2(:,outeriter)=xval;

end

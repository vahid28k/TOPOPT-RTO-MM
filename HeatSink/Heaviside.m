function [rhoH, d_rhoH] = Heaviside(rho_i, beta)

rhoH = 1 - exp(-beta*rho_i) + rho_i*exp(-beta);
%rhoH = rho_i;

if nargout > 1
    d_rhoH = beta*exp(-beta*rho_i) + exp(-beta);
    %d_rhoH = ones(size(rho_i));
end
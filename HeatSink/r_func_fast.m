function [rx, drx] = r_func_fast(x, n, xm, Em, p)

% function computes the coefficients of the elesticity terms used to
% compute the effective young's modulus

% n = number of elements
% Em = number of materials (not including void)

% Initialize arrays

rx = zeros(n, Em+1);
xBlock=reshape(x,n,xm);

% for i = 1:xm
%     xBlock(:, i) = x((i-1)*n+1:i*n);
% end

% RECURSIVE FORMULATION

% drx = zeros(n, Em+1, xm);
% 
% % two materials plus void
% if Em == 2
%     
%     for i = 1:n
%         rx(i, 1) = 1; % void material
%         rx(i, 2) = (xBlock(i, 1)^p)*(xBlock(i, 2)^p); % material 1
%         rx(i, 3) = (xBlock(i, 1)^p)*(1-xBlock(i, 2))^p; % material 2
%     end
%     
% else % Em = 3 -> three materials plus void
%     
%     for i = 1:n
% %         rx(i, 1) = 1; % void material
% %         rx(i, 2) = (xBlock(i, 1)^p)*(xBlock(i, 2)^p)*(xBlock(i, 3)^p); % material 1
% %         rx(i, 3) = (xBlock(i, 1)^p)*(xBlock(i, 2)^p)*(1-xBlock(i, 3))^p; % material 2
% %         rx(i, 4) = (xBlock(i, 1)^p)*(1-xBlock(i, 2))^p; % material 3
%         
%         rx(i, 1) = 1; % void material
%         rx(i, 2) = (xBlock(i, 1)^p)*(xBlock(i, 2)^p)*(xBlock(i, 3)^p); % material 1
%         rx(i, 3) = (xBlock(i, 1)^p)*(xBlock(i, 2)^p)*(1-xBlock(i, 3)^p); % material 2
%         rx(i, 4) = (xBlock(i, 1)^p)*(1-xBlock(i, 2)^p); % material 3
%     end
%     
% end

% SHAPE FUNCTION FORMULATION (Material Interpolation)

drx = zeros(n, (Em+1)* xm);

% used only for 3 material formulation (m=3)


    rx(:,1:(Em+1) ) = [(xBlock(:, 1).^p).*(xBlock(:, 2).^p) ... 
                       (xBlock(:, 1).^p).*((1-xBlock(:, 2)).^p) ... 
                       ((1-xBlock(:, 1)).^p).*(xBlock(:, 2).^p) ... 
                       ((1-xBlock(:, 1)).^p).*((1-xBlock(:, 2)).^p) ];


if nargout > 1  % Compute sensitivities
    
    % RECURSIVE FORMULATION
    
%     if Em == 2
%     
%         for i = 1:n
%             drx(i, 1, :) = zeros(1, 1, 2);
%             
%             drx(i, 2, 1) = p*(xBlock(i, 1)^(p-1))*(xBlock(i, 2)^p);
%             drx(i, 2, 2) = p*(xBlock(i, 2)^(p-1))*(xBlock(i, 1)^p);
%             
%             drx(i, 3, 1) = p*(xBlock(i, 1)^(p-1))*(1-xBlock(i, 2))^p;
%             drx(i, 3, 2) = -(xBlock(i, 1)^p)*p*(1-xBlock(i, 2))^(p-1);
%         end
%         
%     else % Em = 3
%         
%         for i = 1:n
%             drx(i, 1, :) = zeros(1, 1, 3);
%             
%             drx(i, 2, 1) = p*(xBlock(i, 1)^(p-1))*(xBlock(i, 2)^p)*(xBlock(i, 3)^p);
%             drx(i, 2, 2) = p*(xBlock(i, 2)^(p-1))*(xBlock(i, 1)^p)*(xBlock(i, 3)^p);
%             drx(i, 2, 3) = p*(xBlock(i, 3)^(p-1))*(xBlock(i, 1)^p)*(xBlock(i, 2)^p);
%             
%             drx(i, 3, 1) = p*(xBlock(i, 1)^(p-1))*((1-xBlock(i, 3)^p))*(xBlock(i, 2)^p);
%             drx(i, 3, 2) = p*(xBlock(i, 2)^(p-1))*((1-xBlock(i, 3)^p))*(xBlock(i, 1)^p);
%             drx(i, 3, 3) = -(xBlock(i, 1)^p)*(xBlock(i, 2)^p)*(p*xBlock(i, 3)^(p-1));
%             
%             drx(i, 4, 1) = p*(xBlock(i, 1)^(p-1))*(1-xBlock(i, 2)^p);
%             drx(i, 4, 2) = -(xBlock(i, 1)^p)*(p*xBlock(i, 2)^(p-1));
%             drx(i, 4, 3) = 0;           
%         end        
%         
%     end
    
    % SHAPE FUNCTION FORMULATION (Material Interpolation)
    
        drx(:, 1) = p*(xBlock(:, 1).^(p-1)).*(xBlock(:, 2).^p);
        drx(:, 2) = p*(xBlock(:, 2).^(p-1)).*(xBlock(:, 1).^p);

        drx(:, 3) = p*(xBlock(:, 1).^(p-1)).*(1-xBlock(:, 2)).^p;
        drx(:, 4) = -p*(xBlock(:, 1).^p).*(1-xBlock(:, 2)).^(p-1);

        drx(:, 5) = -p*(xBlock(:, 2).^p).*(1-xBlock(:, 1)).^(p-1);
        drx(:, 6) = p*(xBlock(:, 2).^(p-1)).*(1-xBlock(:, 1)).^p;
        
        drx(:, 7) = -p*(1-xBlock(:, 1)).^(p-1).*(1-xBlock(:, 2)).^p;
        drx(:, 8) = -p*(1-xBlock(:, 2)).^(p-1).*(1-xBlock(:, 1)).^p;
            
    
end


        
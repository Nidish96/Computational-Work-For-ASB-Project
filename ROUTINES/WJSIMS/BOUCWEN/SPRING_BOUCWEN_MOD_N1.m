function [Fxyn, varargout] = SPRING_BOUCWEN_MOD_N1(uxyn, pars, opt)
%BOUC-WEN friction model for contact interface
% USAGE:
%  [Fxyn, varargout] = SPRING_BOUCWEN_MOD_N1(uxyn, pars, opt);
% INPUTS:
%   uxyn       : Displacements of contact patch, rows are patches, columns
%                are x,y,z
%   pars       : Parameter set for the contact model [A eta Kp Kn],
%                all are normalized to area. The general Bouc-Wen parameter
%                n is fixed at 1 to allow for analytical solutions.
%                eta = (alpha-beta)
%   opt        : Option matrices for applying parameters to patches and
%                indices on which parameters are log scale. 
% OUTPUTS:
%   Fxyn - Force vectors, rows are partches, columns are x,y,z
%   varargout{1} - Jacobian w.r.t. spatial coordinates
%   varargout{2} - Jacobian w.r.t. parameters

    %convert log parameters to log scale and take derivatives with respect
    %to the conversion. 
    lpars = pars(:);
    lpars(opt.lspci) = 10.^(lpars(opt.lspci));
    dlparsdpars = ones(size(lpars));
    dlparsdpars(opt.lspci) = lpars(opt.lspci).*log(10);    
    
    
    %% Forces    
    Fxyn = zeros(size(uxyn));
    
    % Bouc-Wen - X
    % [A eta Kn]
    A_x = opt.x.T{1}*lpars;
    eta_x = opt.x.T{2}*lpars;
    Kp_x = opt.x.T{3}*lpars;

    ux = abs(uxyn(:, 1));
    sux = sign(uxyn(:, 1));
        
    %calculate the force values
    Fxyn(:, 1) = -sux.*A_x./eta_x .* (exp(-eta_x.*ux) - 1) + Kp_x.*ux.*sux;


    % Bouc-Wen - Y
    % [A eta Kn]
    A_y = opt.y.T{1}*lpars;
    eta_y = opt.y.T{2}*lpars;
    Kp_y = opt.y.T{3}*lpars;
    
    uy = abs(uxyn(:, 2));
    suy = sign(uxyn(:, 2));
        
    %calculate the force values
    Fxyn(:, 2) = -suy.*A_y./eta_y .*(exp(-eta_y.*uy) - 1) + Kp_y.*uy.*suy;
    
    
    % Linear Spring - N
    Kn = opt.n.T{1}*lpars;
    Fxyn(:, 3) = Kn.*uxyn(:,3);
    
    
    %% Jacobians
    if nargout>=2
        
        %this is [partial F_x/partial ux, partial Fy/partial uy, partial
        %Fz/partial uz]
        
        % calculate the set of derivatives: [partial F_x/partial ux, 
        % partial Fy/partial uy, partial Fz/partial uz]
        varargout{1} = [A_x .* exp(-eta_x.*ux) + Kp_x, ...
                        A_y .* exp(-eta_y.*uy) + Kp_y, ...
                        Kn];
    end
    
    if nargout>=3
        varargout{2} = zeros(size(Fxyn,1), length(lpars), 3);
        
        % Pages are X, Y, Z.
        %Columns are [A eta Kn]
        %Rows match displacements still
        
        %Bouc-Wen - X
        dz_dA = -sux./eta_x .* (exp(-eta_x.*ux) - 1);
        
        dz_deta = sux.*A_x./eta_x.^2 .* (exp(-eta_x.*ux) - 1)...
                    -sux.*A_x./eta_x .* (-ux.*exp(-eta_x.*ux));
        
        varargout{2}(:,:,1) = [opt.x.T{1}(:, 1).*dz_dA, ... % A
            opt.x.T{2}(:, 2).*dz_deta, ... %eta
            opt.x.T{3}(:, 3).*uxyn(:, 1),... %Kp
            zeros(size(Fxyn,1), 1)]; % Kn
        
        varargout{2}(:,:,1) = varargout{2}(:,:,1).*dlparsdpars';         
        
        %Bouc-Wen - Y
        dz_dA = -suy./eta_y .* (exp(-eta_y.*uy) - 1);
        
        dz_deta = suy.*A_y./eta_y.^2 .* (exp(-eta_y.*uy) - 1)...
                    -suy.*A_y./eta_y .* (-uy.*exp(-eta_y.*uy));
        
        
        varargout{2}(:,:,2) = [opt.y.T{1}(:, 1).*dz_dA, ... % A
            opt.y.T{2}(:, 2).*dz_deta, ... %eta
            opt.y.T{3}(:, 3).*uxyn(:, 2),... %Kp
            zeros(size(Fxyn,1), 1)]; % Kn
        
        varargout{2}(:,:,2) = varargout{2}(:,:,2).*dlparsdpars';         
        
        
        % Spring - Z
        varargout{2}(:, :, 3) = [zeros(size(uxyn, 1), length(pars)-1) (uxyn(:,3).*opt.n.T{1}(:, 4))].*dlparsdpars';
        
    end

end
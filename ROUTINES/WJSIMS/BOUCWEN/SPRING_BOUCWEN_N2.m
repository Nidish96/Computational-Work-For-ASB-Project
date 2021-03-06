function [Fxyn, varargout] = SPRING_BOUCWEN_N2(uxyn, pars, opt)
%BOUC-WEN friction model for contact interface
% USAGE:
%  [Fxyn, varargout] = SPRING_BOUCWEN_N2(uxyn, pars, opt);
% INPUTS:
%   uxyn       : Displacements of contact patch, rows are patches, columns
%                are x,y,z
%   pars       : Parameter set for the contact model [A eta Kn],
%                all are normalized to area. The general Bouc-Wen parameter
%                n is fixed at 2 to allow for analytical solutions.
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

    u_x = uxyn(:, 1);
        
    %calculate the force values
    Fxyn(:, 1) = sqrt(A_x ./ eta_x) .* tanh(u_x .* sqrt(A_x.*eta_x));


    % Bouc-Wen - Y
    % [A eta Kn]
    A_y = opt.y.T{1}*lpars;
    eta_y = opt.y.T{2}*lpars;
    
    u_y = uxyn(:, 2);
    
    Fxyn(:, 2) = sqrt(A_y ./ eta_y) .* tanh(u_y .* sqrt(A_y.*eta_y));
    
    
    % Linear Spring - N
    Kn = opt.n.T{1}*lpars;
    Fxyn(:, 3) = Kn.*uxyn(:,3);
    
    
    %% Jacobians
    if nargout>=2
        
        %this is [partial F_x/partial ux, partial Fy/partial uy, partial
        %Fz/partial uz]
        
        % calculate the set of derivatives: [partial F_x/partial ux, 
        % partial Fy/partial uy, partial Fz/partial uz]
        varargout{1} = [A_x .* sech(u_x .* sqrt(A_x.*eta_x)).^2, ...
                        A_y .* sech(u_y .* sqrt(A_y.*eta_y)).^2, ...
                        Kn];
    end
    
    if nargout>=3
        varargout{2} = zeros(size(Fxyn,1), length(lpars), 3);
        
        % Pages are X, Y, Z.
        %Columns are [A eta Kn]
        %Rows match displacements still
        
        %Bouc-Wen - X
        dz_dA = 1/2.*sqrt(1./(A_x.*eta_x)).*tanh(u_x .* sqrt(A_x.*eta_x))...
            +u_x./2 .* sech(u_x .* sqrt(A_x.*eta_x)).^2;
        
        dz_deta = -1/2.*sqrt(A_x).*eta_x.^(-3/2).*tanh(u_x .* sqrt(A_x.*eta_x))...
            +A_x ./ eta_x .* u_x./2 .* sech(u_x .* sqrt(A_x.*eta_x)).^2;
        
        
        varargout{2}(:,:,1) = [opt.x.T{1}(:, 1).*dz_dA, ... % A
            opt.x.T{2}(:, 2).*dz_deta, ... %eta
            zeros(size(Fxyn,1), 1)]; % Kn
        
        varargout{2}(:,:,1) = varargout{2}(:,:,1).*dlparsdpars';         
        
        %Bouc-Wen - Y
        dz_dA = 1/2.*sqrt(1./(A_y.*eta_y)).*tanh(u_y .* sqrt(A_y.*eta_y))...
            +u_y./2 .* sech(u_y .* sqrt(A_y.*eta_y)).^2;
        
        dz_deta = -1/2.*sqrt(A_y).*eta_y.^(-3/2).*tanh(u_y .* sqrt(A_y.*eta_y))...
            +A_y ./ eta_y .* u_y./2 .* sech(u_y .* sqrt(A_y.*eta_y)).^2;
        
        
        varargout{2}(:,:,2) = [opt.y.T{1}(:, 1).*dz_dA, ... % A
            opt.y.T{2}(:, 2).*dz_deta, ... %eta
            zeros(size(Fxyn,1), 1)]; % Kn
        
        varargout{2}(:,:,2) = varargout{2}(:,:,2).*dlparsdpars';         
        
        
        % Spring - Z
        varargout{2}(:, :, 3) = [zeros(size(uxyn, 1), length(pars)-1) (uxyn(:,3).*opt.n.T{1}(:, 3))].*dlparsdpars';
        
    end

end
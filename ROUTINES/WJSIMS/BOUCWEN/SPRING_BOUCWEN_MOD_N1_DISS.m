function [Dxyn, varargout] = SPRING_BOUCWEN_MOD_N1_DISS(uxyn, pars, opt)
%JENKINS friction model for contact interface
% USAGE:
%  [Dxyn, varargout] = SPRING_BOUCWEN_MOD_N1_DISS(uxyn, pars, opt);
% INPUTS:
%   uxyn       : Displacements of contact patch, rows are patches, columns
%                are x,y,z
%   pars       : Parameter set for the contact model [A eta Kn],
%                all are normalized to area. The general Bouc-Wen parameter
%                n is fixed at 1 to allow for analytical solutions.
%                eta = (alpha-beta)
%   opt        : Option matrices for applying parameters to patches and
%                indices on which parameters are log scale. 
% OUTPUTS:
%   Dxyn - Dissapation in each direction
%   varargout{1} - Jacobian w.r.t. spatial coordinates 
%   varargout{2} - Jacobian w.r.t. parameters 
    

    %convert log parameters to log scale and take derivatives with respect
    %to the conversion. 
    lpars = pars(:);
    lpars(opt.lspci) = 10.^(lpars(opt.lspci));
    dlparsdpars = ones(size(lpars));
    dlparsdpars(opt.lspci) = lpars(opt.lspci).*log(10); 
    
    %% Forces
    
    %just evaluate the force on the positive side since the goal is the
    %area which should be the same going either direction
    suxyn = sign(uxyn);
    uxyn = abs(uxyn);
    
    %evaluate the force and its derivatives
    %this was in IWAN4, but not used yet here.
    %[Fxyn, dFxyndX, dFxyndp] = SPRING_JENKINS(uxyn, pars, opt);
    
    %dropping third parameter for JENKINS since the derivatives w.r.t
    %parameters are done based on a full loop so secant info is not needed
    [Fxyn, dFxyndX] = SPRING_BOUCWEN_N1(uxyn, pars, opt);
    
    %% Dissipations
    Dxyn = zeros(size(uxyn));
    
    % Bouc-Wen - X
    % [A eta Kn]
    A_x = opt.x.T{1}*lpars;
    eta_x = opt.x.T{2}*lpars;
    
    %displacements in this direction
    ux = uxyn(:,1);
    
    %Applying Masing, dissipation for a full loop.
    Dxyn(:, 1) = (4.*A_x.*ux.*exp(-ux.*eta_x))./eta_x+(8.*A_x.*exp(-ux.*eta_x))./eta_x.^2+(4.*A_x.*ux.*eta_x-8.*A_x)./eta_x.^2;
    
    
    % Bouc-Wen - Y
    % [A eta Kn]
    A_y = opt.y.T{1}*lpars;
    eta_y = opt.y.T{2}*lpars;
    
    %displacements in this direction
    uy = uxyn(:,2);
    
    %Applying Masing, dissipation for a full loop.
    Dxyn(:, 2) = (4.*A_y.*uy.*exp(-uy.*eta_y))./eta_y+(8.*A_y.*exp(-uy.*eta_y))./eta_y.^2+(4.*A_y.*uy.*eta_y-8.*A_y)./eta_y.^2;
    
    
    % Spring - N
    Dxyn(:, 3) = 0;
    
    %% Jacobians
    if nargout>=2
        
        %columns are derivatives in XYZ directions, rows are patches
        %This expression is correct regardless of friction model given that
        %the correct model is called above for Fxyn ect.
        varargout{1} = [4*Fxyn(:, 1:2)-4*dFxyndX(:, 1:2).*uxyn(:, 1:2), ...
            zeros(size(uxyn(:,3)))];
        
        %correct sign output
        varargout{1} = varargout{1}.*suxyn;
    end
    if nargout>=3
        
        %Rows are patches, columns are parameters, pages are XYZ
        %pars are [A eta Kn]
        varargout{2} = zeros(size(Dxyn,1), length(lpars), 3);
        
        % JENKINS 4 - X
        varargout{2}(:, :, 1) = [opt.x.T{1}(:, 1).*( (4.*ux.*exp(-ux.*eta_x))./eta_x+(8.*exp(-ux.*eta_x))./eta_x.^2+(4.*ux.*eta_x-8)./eta_x.^2 ), ... %A
                                 opt.x.T{2}(:, 2).*( -(4.*A_x.*ux.^2.*exp(-ux.*eta_x))./eta_x-(12.*A_x.*ux.*exp(-ux.*eta_x))./eta_x.^2-(16.*A_x.*exp(-ux.*eta_x))./eta_x.^3-(4.*A_x.*ux.*eta_x-16.*A_x)./eta_x.^3 ), ... %eta
                                 zeros(size(uxyn(:, 1:2)))]; %derivative w.r.t. Kp and Kn
        
        varargout{2}(:, :, 1) = varargout{2}(:, :, 1).*dlparsdpars';
        
        %Apply this only if the above calculation is of the forward area
        %not the total area of a full loop. For this one, I am doing a full
        %loop, so this is not needed.
        %varargout{2}(:, :, 1) = -4*dFxyndp(:,:,1).*uxyn(:,1)+8*varargout{2}(:,:,1);

        
        % JENKINS 4 - Y
        varargout{2}(:, :, 2) = [opt.y.T{1}(:, 1).*( (4.*uy.*exp(-uy.*eta_y))./eta_y+(8.*exp(-uy.*eta_y))./eta_y.^2+(4.*uy.*eta_y-8)./eta_y.^2 ), ... %A
                                 opt.y.T{2}(:, 2).*( -(4.*A_y.*uy.^2.*exp(-uy.*eta_y))./eta_y-(12.*A_y.*uy.*exp(-uy.*eta_y))./eta_y.^2-(16.*A_y.*exp(-uy.*eta_y))./eta_y.^3-(4.*A_y.*uy.*eta_y-16.*A_y)./eta_y.^3 ), ... %eta
                                 zeros(size(uxyn(:, 1:2)))]; %derivative w.r.t. Kp and Kn
        
        varargout{2}(:, :, 2) = varargout{2}(:, :, 2).*dlparsdpars';
        
        %Apply this only if the above calculation is of the forward area
        %not the total area of a full loop. For this one, I am doing a full
        %loop, so this is not needed.
        %varargout{2}(:, :, 2) = -4*dFxyndp(:,:,2).*uxyn(:,2)+8*varargout{2}(:,:,2);
        
        % Spring - N: No dissipation
        varargout{2}(:, :, 3) = zeros(size(Dxyn,1), length(lpars));
    end
end
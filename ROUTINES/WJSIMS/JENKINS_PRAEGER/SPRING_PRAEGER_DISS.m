function [Dxyn, varargout] = SPRING_PRAEGER_DISS(uxyn, pars, opt)
%PRAEGER friction model for contact interface
% USAGE:
%  [Dxyn, varargout] = SPRING_PRAEGER_LIN_DISS(uxyn, pars, opt);
% INPUTS:
%   uxyn       : Displacements of contact patch, rows are patches, columns
%                are x,y,z
%   pars       : Parameter set for the contact model [Fs/A; Kp/A; Kn/A]
%   opt        : Option matrices for applying parameters to patches and
%                indices on which parameters are log scale. 
% OUTPUTS:
%   Dxyn - Dissapation in each direction
%   varargout{1} - Jacobian w.r.t. spatial coordinates 
%   varargout{2} - Jacobian w.r.t. parameters 
%
%NOTE: This program leaves Kt as a parameter, but does not pass out related
%derivatives since Kt is fixed when DFUN is set
    

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
%     [Fxyn, dFxyndX] = SPRING_JENKINS(uxyn, pars, opt);
    
    %% Dissipations
    Dxyn = zeros(size(uxyn));
    
    % JENKINS 4 - X
    Fs_x = opt.x.T{1}*lpars;
    
    %displacements in this direction
    ux = uxyn(:,1);
    
    %Applying Masing, dissipation for a full loop. Only applies for
    %macroslip
    Dxyn(:, 1) = 4.*Fs_x.*ux;
    
    
    % JENKINS - Y
    Fs_y = opt.y.T{1}*lpars;
    
    %displacements in this direction
    uy = uxyn(:,2);
    
    %Applying Masing, dissipation for a full loop. Only applies for
    %macroslip
    Dxyn(:, 2) = 4.*Fs_y.*uy;
    
    % Spring - N
    Dxyn(:, 3) = 0;
    
    %% Jacobians
    if nargout>=2
        
        %columns are derivatives in XYZ directions, rows are patches
        %for only in microslip, this expression results in zero. 
        varargout{1} = [4*Fs_x, ...
                        4*Fs_y, ...
            zeros(size(uxyn(:,3)))];
        
        %correct sign output
        varargout{1} = varargout{1}.*suxyn;
    end
    if nargout>=3
        
        %Rows are patches, columns are parameters, pages are XYZ
        %pars are [Fs/A; Kt/A; Kn/A]
        varargout{2} = zeros(size(Dxyn,1), length(lpars), 3);
        
        % JENKINS 4 - X
        varargout{2}(:, :, 1) = [opt.x.T{1}(:, 1).*(4.*ux), ... %Fs
                                 zeros(size(uxyn(:, 1:2)))]; %derivative w.r.t. Kp, Kn
        
        varargout{2}(:, :, 1) = varargout{2}(:, :, 1).*dlparsdpars';
        
        %Apply this only if the above calculation is of the forward area
        %not the total area of a full loop. For Jenkins, I am doing a full
        %loop, so this is not needed.
        %varargout{2}(:, :, 1) = -4*dFxyndp(:,:,1).*uxyn(:,1)+8*varargout{2}(:,:,1);

        
        % JENKINS 4 - Y
        varargout{2}(:, :, 2) = [opt.y.T{1}(:, 1).*(4.*uy), ... %Fs 
                                 zeros(size(uxyn(:, 1:2)))]; %derivative w.r.t. Kp/Kn
        
        varargout{2}(:, :, 2) = varargout{2}(:, :, 2).*dlparsdpars';
        
        %Apply this only if the above calculation is of the forward area
        %not the total area of a full loop. For Jenkins, I am doing a full
        %loop, so this is not needed.
        %varargout{2}(:, :, 2) = -4*dFxyndp(:,:,2).*uxyn(:,2)+8*varargout{2}(:,:,2);
        
        % Spring - N: No dissipation
        varargout{2}(:, :, 3) = zeros(size(Dxyn,1), length(lpars));
    end
end
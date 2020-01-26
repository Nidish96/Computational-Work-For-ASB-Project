function [Fxyn, varargout] = SPRING_PRAEGER(uxyn, pars, opt)
%PRAEGER friction model for contact interface
% USAGE:
%  [Fxyn, varargout] = SPRING_PRAEGER_LIN(uxyn, pars, opt);
% INPUTS:
%   uxyn       : Displacements of contact patch, rows are patches, columns
%                are x,y,z
%   pars       : Parameter set for the contact model [Fs/A; Kp/A; Kn/A]
%   opt        : Option matrices for applying parameters to patches and
%                indices on which parameters are log scale. 
% OUTPUTS:
%   Fxyn - Force vectors, rows are partches, columns are x,y,z
%   varargout{1} - Jacobian w.r.t. spatial coordinates
%   varargout{2} - Jacobian w.r.t. parameters
%
%NOTE: This program leaves Kt as a parameter, but does not pass out related
%derivatives since Kt is fixed when CFUN is set

    %convert log parameters to log scale and take derivatives with respect
    %to the conversion. 
    lpars = pars(:);
    lpars(opt.lspci) = 10.^(lpars(opt.lspci));
    dlparsdpars = ones(size(lpars));
    dlparsdpars(opt.lspci) = lpars(opt.lspci).*log(10);    
    
    
    %% Forces    
    Fxyn = zeros(size(uxyn));
    
    % Jenkins - X
    Fs = opt.x.T{1}*lpars;
    Kp = opt.x.T{2}*lpars;
    
    sux = sign(real(uxyn(:,1)));
    ux = sux.*uxyn(:,1);
    Fxyn(:,1) = Fs.*sux + Kp.*ux.*sux;
    
    % Jenkins - Y
    fs = opt.y.T{1}*lpars;
    kp = opt.y.T{2}*lpars;
    
    suy = sign(real(uxyn(:,2)));
    uy = suy.*uxyn(:,2);
	Fxyn(:,2) = fs.*suy + kp.*uy.*suy;
    
    % Linear Spring - N
    Kn = opt.n.T{1}*lpars;
    Fxyn(:, 3) = Kn.*uxyn(:,3);
    
    %% Jacobians
    if nargout>=2
        
        %this is [partial F_x/partial ux, partial Fy/partial uy, partial
        %Fz/partial uz]
        
        varargout{1} = [Kp, kp, Kn];
    end
    if nargout>=3
        varargout{2} = zeros(size(Fxyn,1), length(lpars), 3);
        
        %Jenkins - X
        varargout{2}(:,:,1) = [opt.x.T{1}(:, 1).*sux, ... %Fs
            opt.x.T{2}(:, 2).*uxyn(:, 1), ... %Kp
            zeros(size(Fxyn,1), 1)];
        
        varargout{2}(:,:,1) = varargout{2}(:,:,1).*dlparsdpars';  % Check        
        
        %Jenkins - Y
        varargout{2}(:,:,2) = [opt.y.T{1}(:, 1).*suy, ... %fs
            opt.y.T{2}(:, 2).*uxyn(:, 2), ... %Kp
            zeros(size(Fxyn,1), 1)];
        varargout{2}(:,:,2) = varargout{2}(:,:,2).*dlparsdpars';  % Check        
        
        % Spring - Z
        varargout{2}(:, :, 3) = [zeros(size(uxyn, 1), length(pars)-1) (uxyn(:,3).*opt.n.T{1}(:, 3))].*dlparsdpars';
    end
end
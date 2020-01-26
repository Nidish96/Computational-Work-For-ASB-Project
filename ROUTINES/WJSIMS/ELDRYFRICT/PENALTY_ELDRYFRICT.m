function [Fxyn, varargout] = PENALTY_ELDRYFRICT(uxyn, pars, opt)
%JENKINS friction model for contact interface
% USAGE:
%  [Fxyn, varargout] = PENALTY_ELDRYFRICT(uxyn, pars, opt);
% INPUTS:
%   uxyn       : Displacements of contact patch, rows are patches, columns
%                are x,y,z
%   pars       : Parameter set for the contact model [mu; Kt/A; Kn/A]
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
    
    % Penalty Linear Spring - N
    Kn = opt.n.T{1}*lpars;
    Fxyn(:, 3) = max(Kn.*uxyn(:,3),0);
    coni = (Fxyn(:,3)~=0);
    
    % Jenkins - X
    Mu = opt.x.T{1}*lpars;
    Fs = Mu.*Fxyn(:,3);
    Kt = opt.x.T{2}*lpars;
    
    sux = sign(real(uxyn(:,1)));
    ux = sux.*uxyn(:,1);
    Fxyn(:,1) = (Kt.*ux.*(Kt.*ux < Fs) + Fs.*(Kt.*ux >= Fs)).*sux.*coni;
    
    % Jenkins - Y
    mu = opt.y.T{1}*lpars;
    fs = mu.*Fxyn(:,3);
    kt = opt.y.T{2}*lpars;
    
    suy = sign(real(uxyn(:,2)));
    uy = suy.*uxyn(:,2);
	Fxyn(:,2) = (kt.*uy.*(kt.*uy < fs) + fs.*(kt.*uy >= fs)).*suy.*coni;
    
    %% Jacobians
    if nargout>=2
        
        %this is [partial F_x/partial ux, partial Fy/partial uy, partial
        %Fz/partial uz]
        varargout{1} = zeros(size(Fxyn,1),3,3);
        % X Force
        varargout{1}(:,1,1) = Kt.*(Kt.*ux < Fs).*coni;
        varargout{1}(:,3,1) = Mu.*Kn.*(Kt.*ux >= Fs).*coni.*sux;
        % Y Force
        varargout{1}(:,2,2) = kt.*(kt.*uy < fs).*coni;
        varargout{1}(:,3,2) = mu.*Kn.*(kt.*uy >= fs).*coni.*suy;
        % N Force
        varargout{1}(:,3,3) = Kn.*coni;
        
%         varargout{1} = [Kt.*(Kt.*ux < Fs), kt.*(kt.*uy < fs), Kn].*coni;
    end
    if nargout>=3
        varargout{2} = zeros(size(Fxyn,1), length(lpars), 3);
        
        %Jenkins - X
        varargout{2}(:,:,1) = [opt.x.T{1}(:, 1).*Fxyn(:,3).*sux.*(Kt.*ux >= Fs), ...
            opt.x.T{2}(:, 2).*sux.*ux.*(Kt.*ux < Fs), ...
            opt.n.T{1}(:, 3).*(Mu.*uxyn(:,3)).*sux.*(Kt.*ux >= Fs) %Justin's
            ];
        %             2*opt.x.T{1}(:,1).*((opt.n.T{1}(:,3)*lpars(3)).*uxyn(:,3)).*sux.*(Kt.*ux >= Fs) %Nidish's last row. 
        
        varargout{2}(:,:,1) = varargout{2}(:,:,1).*coni.*dlparsdpars';  % Check        
        
        %Jenkins - Y
        varargout{2}(:,:,2) = [opt.y.T{1}(:, 1).*Fxyn(:,3).*suy.*(kt.*uy >= fs), ...
            opt.y.T{2}(:, 2).*suy.*uy.*(kt.*uy < fs),...
            opt.n.T{1}(:, 3).*(mu.*uxyn(:,3)).*suy.*(kt.*uy >= fs)];
%         2*opt.y.T{1}(:,1).*((opt.n.T{1}(:,3)*lpars(3)).*uxyn(:,3)).*suy.*(kt.*uy >= fs) %Nidish's last row. 

        varargout{2}(:,:,2) = varargout{2}(:,:,2).*coni.*dlparsdpars';  % Check        
        
        % Spring - Z
        varargout{2}(:, :, 3) = [zeros(size(uxyn, 1), length(pars)-1) (uxyn(:,3).*opt.n.T{1}(:, 3)).*coni].*dlparsdpars';
    end
end
function [Fxyn, varargout] = SPRING_IWAN5(uxyn, pars, opt)
%IWAN-5 friction model for contact interface
% USAGE:
%  [Fxyn, varargout] = SPRING_IWAN5(uxyn, pars, opt);
% INPUTS:
%   uxyn       : Displacements of contact patch, rows are patches, columns
%                are x,y,z
%   pars       : Parameter set for the contact model [Fs kt Chi bt theta],
%                normalized to area? theta is not normalized to area.  
%   opt        : Option matrices for applying parameters to patches and
%                indices on which parameters are log scale. 
% OUTPUTS:
%   Fxyn - Force vectors, rows are partches, columns are x,y,z
%   varargout{1} - Jacobian w.r.t. spatial coordinates
%   varargout{2} - Jacobian w.r.t. parameters pages are for x,y,z
%               directions
%disp(mat2str(pars));

    lpars = pars(:);
    lpars(opt.lspci) = 10.^(lpars(opt.lspci));
    dlparsdpars = ones(size(lpars));
    dlparsdpars(opt.lspci) = lpars(opt.lspci).*log(10);    
    
    %% Forces    
    Fxyn = zeros(size(uxyn));
    
    % Iwan 4 - X
    Fs_x = opt.x.T{1}*lpars;
    Kt_x = opt.x.T{2}*lpars;
    Chi_x = opt.x.T{3}*lpars;
    Bt_x = opt.x.T{4}*lpars;
    Theta_x = opt.x.T{5}*lpars;
    
    PhiMx_x = Fs_x.*(1+Bt_x)./(Kt_x.*(Bt_x+(Chi_x+1)./(Chi_x+2)));
    CoefE_x = (Kt_x.*(Kt_x.*(Bt_x+(Chi_x+1)./(Chi_x+2))./(Fs_x.*(1+Bt_x))).^(1+Chi_x))./((1+Bt_x).*(Chi_x+2)).*(Theta_x.*(Chi_x + 1) - Chi_x - 2);    
    
    sux = sign(real(uxyn(:,1)));
    ux = sux.*uxyn(:,1);
	Fxyn(:,1) = ((Kt_x.*ux + CoefE_x.*ux.^(2+Chi_x)).*(ux<PhiMx_x) + (Fs_x.*Theta_x).*(ux>=PhiMx_x)).*sux;
    
    % Iwan 4 - Y
    Fs_y = opt.y.T{1}*lpars;
    Kt_y = opt.y.T{2}*lpars;
    Chi_y = opt.y.T{3}*lpars;
    Bt_y = opt.y.T{4}*lpars;
    Theta_y = opt.y.T{5}*lpars;
    
    PhiMx_y = Fs_y.*(1+Bt_y)./(Kt_y.*(Bt_y+(Chi_y+1)./(Chi_y+2)));
    CoefE_y = (Kt_y.*(Kt_y.*(Bt_y+(Chi_y+1)./(Chi_y+2))./(Fs_y.*(1+Bt_y))).^(1+Chi_y))./((1+Bt_y).*(Chi_y+2)).*(Theta_y.*(Chi_y + 1) - Chi_y - 2);
   
    suy = sign(real(uxyn(:,2)));
    uy = suy.*uxyn(:,2);
	Fxyn(:,2) = ((Kt_y.*uy + CoefE_y.*uy.^(2+Chi_y)).*(uy<PhiMx_y) + (Fs_y.*Theta_y).*(uy>=PhiMx_y)).*suy;
    
    % Linear Spring - N
    Kn = opt.n.T{1}*lpars;
    Fxyn(:, 3) = Kn.*uxyn(:,3);
    %% Jacobians
    if nargout>=2
        varargout{1} = [(Kt_x+CoefE_x.*(2+Chi_x).*ux.^(1+Chi_x)).*(ux<PhiMx_x), (Kt_y+CoefE_y.*(2+Chi_y).*uy.^(1+Chi_y)).*(uy<PhiMx_y), Kn];
    end
    if nargout>=3
        varargout{2} = zeros(size(Fxyn,1), length(lpars), 3);
        
        % Iwan 4 - X
%         temp_num = ((Bt_x+1).*Chi_x.^2+((-Bt_x-1).*Theta_x+4.*Bt_x+3).*Chi_x+(-2.*Bt_x-1).*Theta_x+4.*Bt_x+2)...
%             .*log(((Kt_x.*Bt_x+Kt_x).*Chi_x+2.*Kt_x.*Bt_x+Kt_x)./((Fs_x.*Bt_x+Fs_x).*Chi_x+2.*Fs_x.*Bt_x+2.*Fs_x))...
%             +Chi_x+Bt_x.*Theta_x+1;
%         temp_den = (Bt_x+1).*Chi_x.^2+((-Bt_x-1).*Theta_x+4.*Bt_x+3).*Chi_x+(-2.*Bt_x-1).*Theta_x+4.*Bt_x+2;
        
        temp_num = (((Bt_x+1).*Theta_x-Bt_x-1).*Chi_x.^2+((3.*Bt_x+2).*Theta_x-4.*Bt_x-3).*Chi_x+(2.*Bt_x+1).*Theta_x-4.*Bt_x-2)...
            .*log(((Kt_x.*Bt_x+Kt_x).*Chi_x+2.*Kt_x.*Bt_x+Kt_x)./((Fs_x.*Bt_x+Fs_x).*Chi_x+2.*Fs_x.*Bt_x+2.*Fs_x))...
            +(Theta_x-1).*Chi_x+(Bt_x+1).*Theta_x-1;
        temp_den = ((Bt_x+1).*Theta_x-Bt_x-1).*Chi_x.^2+((3.*Bt_x+2).*Theta_x-4.*Bt_x-3).*Chi_x...
            +(2.*Bt_x+1).*Theta_x-4.*Bt_x-2;
        
        dCdp = CoefE_x.*[-(1+Chi_x)./Fs_x, ... %Fs [Fs kt Chi bt theta]
                        (2+Chi_x)./Kt_x, ... %Kt
                        temp_num./temp_den, ... % chi - not implemented yet
                        -(Bt_x.*(Chi_x+2))./((Bt_x+1).*(Bt_x.*Chi_x+Chi_x+2.*Bt_x+1)), ... % beta
                        (Chi_x + 1) ./ (Theta_x .* (Chi_x + 1) - Chi_x - 2) % theta
                        ]; 
                    
        tmp = [dCdp(:,1).*ux.^(2+Chi_x), ...
                ux+dCdp(:,2).*ux.^(2+Chi_x), ...
                (dCdp(:,3)+CoefE_x.*log(ux)).*ux.^(Chi_x+2), ...
                dCdp(:,4).*ux.^(2+Chi_x), ...
                dCdp(:, 5).*ux.^(2+Chi_x)
                ].*(ux<PhiMx_x) ...
             + [Theta_x zeros(length(Theta_x), 3) Fs_x].*(ux>=PhiMx_x);
        tmp = sux.*tmp;
        for i=1:size(tmp,2)
            varargout{2}(:, :, 1) = varargout{2}(:, :, 1) + diag(tmp(:,i))*opt.x.T{i};
        end
        varargout{2}(:,:,1) = varargout{2}(:,:,1).*dlparsdpars';  % Check
        
        % Iwan 4 - Y
%         t1 = log(Kt_y.*(Bt_y+(Chi_y+1)./(Chi_y+2))./(Fs_y.*(Bt_y+1)));
%         t2 = Bt_y./((Bt_y+1).*Chi_y+2*Bt_y+1);
%         dCdp = coefE.*[-(1+Chi_y)./Fs_y, ... %Fs [Fs kt Chi bt theta]
%                         (2+Chi_y)./Kt_y, ... %Kt
%                         t1-t2, ... % Chi
%                         -Bt_y./((Bt_y+1).*(Bt_y+(1+Chi_y)./(2+Chi_y)))% bt
%                         ];
%             
%         tmp = [-dCdp(:,1).*uy.^(2+Chi_y), ...
%                 uy-dCdp(:,2).*uy.^(2+Chi_y), ...
%               -(dCdp(:,3)+coefE.*log(uy)).*uy.^(Chi_y+2), ...
%                -dCdp(:,4).*uy.^(2+Chi_y)].*(uy<phiMx) + ...
%              [1 0 0 0].*(uy>=phiMx);
%         tmp = suy.*tmp;
%         temp_num = ((Bt_y+1).*Chi_y.^2+((-Bt_y-1).*Theta_y+4.*Bt_y+3).*Chi_y+(-2.*Bt_y-1).*Theta_y+4.*Bt_y+2)...
%             .*log(((Kt_y.*Bt_y+Kt_y).*Chi_y+2.*Kt_y.*Bt_y+Kt_y)./((Fs_y.*Bt_y+Fs_y).*Chi_y+2.*Fs_y.*Bt_y+2.*Fs_y))...
%             +Chi_y+Bt_y.*Theta_y+1;
%         temp_den = (Bt_y+1).*Chi_y.^2+((-Bt_y-1).*Theta_y+4.*Bt_y+3).*Chi_y+(-2.*Bt_y-1).*Theta_y+4.*Bt_y+2;
        
        temp_num = (((Bt_y+1).*Theta_y-Bt_y-1).*Chi_y.^2+((3.*Bt_y+2).*Theta_y-4.*Bt_y-3).*Chi_y+(2.*Bt_y+1).*Theta_y-4.*Bt_y-2)...
            .*log(((Kt_y.*Bt_y+Kt_y).*Chi_y+2.*Kt_y.*Bt_y+Kt_y)./((Fs_y.*Bt_y+Fs_y).*Chi_y+2.*Fs_y.*Bt_y+2.*Fs_y))...
            +(Theta_y-1).*Chi_y+(Bt_y+1).*Theta_y-1;
        temp_den = ((Bt_y+1).*Theta_y-Bt_y-1).*Chi_y.^2+((3.*Bt_y+2).*Theta_y-4.*Bt_y-3).*Chi_y...
            +(2.*Bt_y+1).*Theta_y-4.*Bt_y-2;
        
        dCdp = CoefE_y.*[-(1+Chi_y)./Fs_y, ... %Fs [Fs kt Chi bt theta]
                        (2+Chi_y)./Kt_y, ... %Kt
                        temp_num./temp_den, ... % chi - not implemented yet
                        -(Bt_y.*(Chi_y+2))./((Bt_y+1).*(Bt_y.*Chi_y+Chi_y+2.*Bt_y+1)), ... % beta
                        (Chi_y + 1) ./ (Theta_y .* (Chi_y + 1) - Chi_y - 2) % theta
                        ]; 
                    
        tmp = [dCdp(:,1).*uy.^(2+Chi_y), ...
                uy+dCdp(:,2).*uy.^(2+Chi_y), ...
                (dCdp(:,3)+CoefE_y.*log(uy)).*uy.^(Chi_y+2), ...
                dCdp(:,4).*uy.^(2+Chi_y), ...
                dCdp(:, 5).*uy.^(2+Chi_y)
                ].*(uy<PhiMx_y) ...
             + [Theta_y zeros(length(Theta_y), 3) Fs_y].*(uy>=PhiMx_y);
         
        tmp = suy.*tmp;
        for i=1:size(tmp,2)
            varargout{2}(:, :, 2) = varargout{2}(:, :, 2) + diag(tmp(:,i))*opt.y.T{i};
        end
        varargout{2}(:,:,2) = varargout{2}(:,:,2).*dlparsdpars';  % Check
        
        % Penalty - Z
        tmp = uxyn(:,3);
        varargout{2}(:, :, 3) = (diag(tmp)*opt.n.T{1}).*dlparsdpars';
    end
end
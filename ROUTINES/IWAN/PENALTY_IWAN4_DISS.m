function [Dxyn, varargout] = PENALTY_IWAN4_DISS(uxyn, pars, opt)
    lpars = pars;
    lpars(opt.lspci) = 10.^(lpars(opt.lspci));
    dlparsdpars = ones(size(pars));
    dlparsdpars(opt.lspci) = lpars(opt.lspci).*log(10); 
    
    %% Forces
    suxyn = sign(uxyn);
    uxyn = abs(uxyn);
    [Fxyn, dFxyndX, dFxyndp] = PENALTY_IWAN4(uxyn, pars, opt);
    
    %% Dissipations
    Dxyn = zeros(size(uxyn));
    % Iwan 4 - X
    Fs = opt.x.T{1}*lpars;
    Kt = opt.x.T{2}*lpars;
    Chi = opt.x.T{3}*lpars;
    Bt = opt.x.T{4}*lpars;
    
    PhiMx = Fs.*(1+Bt)./(Kt.*(Bt+(Chi+1)./(Chi+2)));
    CoefE = (Kt.*(Kt.*(Bt+(Chi+1)./(Chi+2))./(Fs.*(1+Bt))).^(1+Chi))./((1+Bt).*(Chi+2));
    
    ux = uxyn(:,1);
    Dxyn(:, 1) = (Kt/2.*ux.^2-CoefE.*ux.^(Chi+3)./(Chi+3)).*(ux<PhiMx) + ...
        (Kt/2.*PhiMx.^2-CoefE.*PhiMx.^(Chi+3)./(Chi+3)+Fs.*(ux-PhiMx)).*(ux>=PhiMx);
    Dxyn(:, 1) = -4*Fxyn(:,1).*ux+8*Dxyn(:,1);
    
    Dxyn(:, 1) = (4*CoefE.*ux.^(Chi+3).*(Chi+1)./(Chi+3)).*(ux<PhiMx) + ...
        (4*CoefE.*PhiMx.^(Chi+3).*(Chi+1)./(Chi+3)+4*Fs.*PhiMx.*ux).*(ux>=PhiMx);
    
    % Iwan 4 - Y
    fs = opt.y.T{1}*lpars;
    kt = opt.y.T{2}*lpars;
    chi = opt.y.T{3}*lpars;
    bt = opt.y.T{4}*lpars;
    
    phiMx = fs.*(1+bt)./(kt.*(bt+(chi+1)./(chi+2)));
    coefE = (kt.*(kt.*(bt+(chi+1)./(chi+2))./(fs.*(1+bt))).^(1+chi))./((1+bt).*(chi+2));    
    
    uy = uxyn(:,2);
    Dxyn(:, 2) = (kt/2.*uy.^2-coefE.*uy.^(chi+3)./(chi+3)).*(uy<phiMx) + ...
        (kt/2.*phiMx.^2-coefE.*phiMx.^(chi+3)./(chi+3)+fs.*(uy-phiMx)).*(uy>=phiMx);
    Dxyn(:, 2) = -4*Fxyn(:,2).*uy+8*Dxyn(:,2);
    
    Dxyn(:, 2) = (4*coefE.*uy.^(chi+3).*(chi+1)./(chi+3)).*(uy<phiMx) + ...
        (4*coefE.*phiMx.^(chi+3).*(chi+1)./(chi+3)+4*fs.*phiMx.*uy).*(uy>=phiMx);
    
    % Penalty - N
    Dxyn(:, 3) = 0;
    
    %% Jacobians
    if nargout>=2
        varargout{1} = [(Kt.*ux-CoefE.*ux.^(Chi+2)).*(ux<PhiMx)+(Fs).*(ux>=PhiMx), ...
            (kt.*uy-coefE.*uy.^(chi+2)).*(uy<phiMx)+(fs).*(uy>=phiMx), ...
            zeros(size(uxyn(:,3)))];
        varargout{1}(:, 1:2) = -4*dFxyndX(:, 1:2).*uxyn(:, 1:2)-4*Fxyn(:, 1:2)+8*varargout{1}(:, 1:2);
        varargout{1} = varargout{1}.*suxyn;
    end
    if nargout>=3
        varargout{2} = zeros(size(Dxyn,1), length(lpars), 3);
        
        % Iwan 4 - X
        t1 = log(Kt.*(Bt+(Chi+1)./(Chi+2))./(Fs.*(Bt+1)));
        t2 = Bt./((Bt+1).*Chi+2*Bt+1);
        dCdp = CoefE.*[-(1+Chi)./Fs, ...
                        (2+Chi)./Kt, ...
                        t1-t2, ...
                        (-Bt./((Bt+1).*(Bt+(1+Chi)./(2+Chi))))];
        dPhidp = PhiMx.*[1./Fs, ...
            -1./Kt, ...
            -1./((Chi+2).^2.*(Bt+(Chi+1)./(Chi+2))), ...
            -1./((Bt+1).*(Chi+2).*(Bt+(Chi+1)./(Chi+2)))];
        tmp = ([-dCdp(:,1).*(ux.^(Chi+3)./(Chi+3)), ...
            ux.^2/2-dCdp(:,2).*(ux.^(Chi+3)./(Chi+3)), ...
            -dCdp(:,3).*(ux.^(Chi+3)./(Chi+3))+CoefE.*(ux.^(Chi+3)./(Chi+3)).*(1./(Chi+3)-log(ux)), ...
            -dCdp(:,4).*(ux.^(Chi+3)./(Chi+3))]).*(ux<PhiMx) + ...
            ([(Kt.*PhiMx-CoefE.*PhiMx.^(Chi+2)-Fs).*dPhidp(:,1) + (-(PhiMx.^(Chi+3)./(Chi+3))+(ux-PhiMx)).*dCdp(:,1), ...
            PhiMx.^2/2-(Kt.*PhiMx-CoefE.*PhiMx.^(Chi+2)-Fs).*dPhidp(:,2) + (-(PhiMx.^(Chi+3)./(Chi+3))+(ux-PhiMx)).*dCdp(:,2), ...
            (Kt.*PhiMx-CoefE.*PhiMx.^(Chi+2)-Fs).*dPhidp(:,3) + (-(PhiMx.^(Chi+3)./(Chi+3))+(ux-PhiMx)).*dCdp(:,3) - CoefE.*PhiMx.^(Chi+3).*log(PhiMx)./(Chi+3), ...
            (Kt.*PhiMx-CoefE.*PhiMx.^(Chi+2)-Fs).*dPhidp(:,4) + (-(PhiMx.^(Chi+3)./(Chi+3))+(ux-PhiMx)).*dCdp(:,4)]).*(ux>=PhiMx);
        for i=1:size(tmp,2)
            varargout{2}(:, :, 1) = varargout{2}(:, :, 1) + diag(tmp(:,i))*opt.x.T{i};
        end
        varargout{2}(:, :, 1) = varargout{2}(:, :, 1).*dlparsdpars';
        varargout{2}(:, :, 1) = -4*dFxyndp(:,:,1).*uxyn(:,1)+8*varargout{2}(:,:,1);

        % Iwan 4 - Y
        t1 = log(kt.*(bt+(chi+1)./(chi+2))./(fs.*(bt+1)));
        t2 = bt./((bt+1).*chi+2*bt+1);
        dcdp = coefE.*[-(1+chi)./fs, ...
                        (2+chi)./kt, ...
                        t1-t2, ...
                        (-bt./((bt+1).*(bt+(1+chi)./(2+chi))))];
        dPhidp = phiMx.*[1./fs, ...
            -1./kt, ...
            -1./((chi+2).^2.*(bt+(chi+1)./(chi+2))), ...
            -1./((bt+1).*(chi+2).*(bt+(chi+1)./(chi+2)))];
        tmp = ([-dcdp(:,1).*(uy.^(chi+3)./(chi+3)), ...
            uy.^2/2-dcdp(:,2).*(uy.^(chi+3)./(chi+3)), ...
            -dcdp(:,3).*(uy.^(chi+3)./(chi+3))+coefE.*(uy.^(chi+3)./(chi+3)).*(1./(chi+3)-log(uy)), ...
            -dcdp(:,4).*(uy.^(chi+3)./(chi+3))]).*(uy<phiMx) + ...
            ([(kt.*phiMx-coefE.*phiMx.^(chi+2)-fs).*dPhidp(:,1) + (-(phiMx.^(chi+3)./(chi+3))+(uy-phiMx)).*dcdp(:,1), ...
            phiMx.^2/2-(kt.*phiMx-coefE.*phiMx.^(chi+2)-fs).*dPhidp(:,2) + (-(phiMx.^(chi+3)./(chi+3))+(uy-phiMx)).*dcdp(:,2), ...
            (kt.*phiMx-coefE.*phiMx.^(chi+2)-fs).*dPhidp(:,3) + (-(phiMx.^(chi+3)./(chi+3))+(uy-phiMx)).*dcdp(:,3) - coefE.*phiMx.^(chi+3).*log(phiMx)./(chi+3), ...
            (kt.*phiMx-coefE.*phiMx.^(chi+2)-fs).*dPhidp(:,4) + (-(phiMx.^(chi+3)./(chi+3))+(uy-phiMx)).*dcdp(:,4)]).*(uy>=phiMx);
        for i=1:size(tmp,2)
            varargout{2}(:, :, 2) = varargout{2}(:, :, 2) + diag(tmp(:,i))*opt.y.T{i};
        end
        varargout{2}(:, :, 2) = varargout{2}(:, :, 2).*dlparsdpars';
        varargout{2}(:, :, 2) = -4*dFxyndp(:,:,2).*uxyn(:,2)+8*varargout{2}(:,:,2);
        
        % Penalty - N: No dissipation
        varargout{2}(:, :, 3) = zeros(size(Dxyn,1), length(lpars));
    end
end
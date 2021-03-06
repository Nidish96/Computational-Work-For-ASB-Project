function [Dxyn, varargout] = SPRING_IWAN4_DISS(uxyn, pars, opt)
%IWAN 4 friction model for contact interface - I STARTED MESSING WITH THIS
%FUNCTION, NO GUARENTEES THAT IT IS CORRECT. 
% USAGE:
%  [Dxyn, varargout] = SPRING_BOUCWEN_DISS(uxyn, pars, opt);
% INPUTS:
%   uxyn       : Displacements of contact patch, rows are patches, columns
%                are x,y,z
%   pars       : Parameter set for the contact model [A eta n Kn],
%                all except for n are normalized to area. 
%   opt        : Option matrices for applying parameters to patches and
%                indices on which parameters are log scale. 
% OUTPUTS:
%   Dxyn - Dissapation in each direction
%   varargout{1} - Jacobian w.r.t. spatial coordinates
%   varargout{2} - Jacobian w.r.t. parameters

    lpars = pars(:);
    lpars(opt.lspci) = 10.^(lpars(opt.lspci));
    dlparsdpars = ones(size(pars));
    dlparsdpars(opt.lspci) = lpars(opt.lspci).*log(10); 
    
    %% Forces
    suxyn = sign(uxyn);
    uxyn = abs(uxyn);
    [Fxyn, dFxyndX, dFxyndp] = SPRING_IWAN4(uxyn, pars, opt);
    
    %% Dissipations
    Dxyn = zeros(size(uxyn));
    
    % Iwan 4 - X
    Fs_x = opt.x.T{1}*lpars;
    Kt_x = opt.x.T{2}*lpars;
    Chi_x = opt.x.T{3}*lpars;
    Bt_x = opt.x.T{4}*lpars;
    
    PhiMx_x = Fs_x.*(1+Bt_x)./(Kt_x.*(Bt_x+(Chi_x+1)./(Chi_x+2)));
    
    ux = uxyn(:,1);
    
    %Stick - potentially wrong (had missing multiply in wxMaxima, may have
    %misinterpreted. 
    D1 = (4.*ux.^(Chi_x+3).*(Fs_x.*(Bt_x+1)./(Kt_x.*((Chi_x+1)./(Chi_x+2)+Bt_x))).^(-Chi_x-2).*Fs_x.*(Chi_x+1))...
        ./(((Chi_x+1)./(Chi_x+2)+Bt_x).*(Chi_x+2).*(Chi_x+3));
    
    %Slip
    D2 = (4.*(ux./(Chi_x+2)-Fs_x.*(Bt_x+1)./(Kt_x.*(Chi_x+3).*((Chi_x+1)./(Chi_x+2)+Bt_x))).*Fs_x.*(Chi_x+1))...
        ./((Chi_x+1)./(Chi_x+2)+Bt_x)...
        +(4.*Fs_x.*Bt_x.*(ux-Fs_x.*(Bt_x+1)./(Kt_x.*((Chi_x+1)./(Chi_x+2)+Bt_x))))...
        ./((Chi_x+1)./(Chi_x+2)+Bt_x);
    
    Dxyn(:, 1) = D1.*(ux<PhiMx_x) + D2.*(ux>=PhiMx_x);
    
    % Iwan 4 - Y
    Fs_y = opt.y.T{1}*lpars;
    Kt_y = opt.y.T{2}*lpars;
    Chi_y = opt.y.T{3}*lpars;
    Bt_y = opt.y.T{4}*lpars;
    
    PhiMx_y = Fs_y.*(1+Bt_y)./(Kt_y.*(Bt_y+(Chi_y+1)./(Chi_y+2)));   
    
    uy = uxyn(:,2);
    D1 = (4.*uy.^(Chi_y+3).*(Fs_y.*(Bt_y+1)./(Kt_y.*((Chi_y+1)./(Chi_y+2)+Bt_y))).^(-Chi_y-2).*Fs_y.*(Chi_y+1))...
        ./(((Chi_y+1)./(Chi_y+2)+Bt_y).*(Chi_y+2).*(Chi_y+3));
    
    %Slip
    D2 = (4.*(uy./(Chi_y+2)-Fs_y.*(Bt_y+1)./(Kt_y.*(Chi_y+3).*((Chi_y+1)./(Chi_y+2)+Bt_y))).*Fs_y.*(Chi_y+1))...
        ./((Chi_y+1)./(Chi_y+2)+Bt_y)...
        +(4.*Fs_y.*Bt_y.*(uy-Fs_y.*(Bt_y+1)./(Kt_y.*((Chi_y+1)./(Chi_y+2)+Bt_y))))...
        ./((Chi_y+1)./(Chi_y+2)+Bt_y);
    
    Dxyn(:, 2) = D1.*(uy<PhiMx_y) + D2.*(uy>=PhiMx_y);
    
    % Penalty - N
    Dxyn(:, 3) = 0;
    
    %% Jacobians
    if nargout>=2
%         varargout{1} = [(Kt_x.*ux-CoefE_x.*ux.^(Chi_x+2)).*(ux<PhiMx_x)+(Fs_x).*(ux>=PhiMx_x), ...
%             (kt.*uy-coefE.*uy.^(chi+2)).*(uy<phiMx)+(fs).*(uy>=phiMx), ...
%             zeros(size(uxyn(:,3)))];
%         varargout{1}(:, 1:2) = -4*dFxyndX(:, 1:2).*uxyn(:, 1:2)-4*Fxyn(:, 1:2)+8*varargout{1}(:, 1:2);
%         varargout{1} = varargout{1}.*suxyn;

        %columns are derivatives in XYZ directions, rows are patches
        %for only in microslip, this expression results in zero. 
        varargout{1} = [4*Fxyn(:, 1:2)-4*dFxyndX(:, 1:2).*uxyn(:, 1:2), ...
            zeros(size(uxyn(:,3)))];
        
        %correct sign output
        varargout{1} = varargout{1}.*suxyn;
        
        %Note that this does have a bit of numerical error just because the
        %two components are very close in magnitude and then subtracted. 
        
    end
    if nargout>=3
        varargout{2} = zeros(size(Dxyn,1), length(lpars), 3);
        
        % Iwan 4 - X [Fs Kt Chi Bt]
    
        dD1dFs = -(4.*Kt_x.^2.*ux.^(Chi_x+3).*(Chi_x+1).^2.*(Bt_x.*Chi_x+Chi_x+2.*Bt_x+1))...
            ./(Fs_x.^2.*(Bt_x+1).^2.*(Chi_x+2).^2.*(Chi_x+3).*((Fs_x.*(Bt_x+1).*(Chi_x+2))...
            ./(Kt_x.*((Bt_x+1).*Chi_x+2.*Bt_x+1))).^Chi_x);
        
        dD1dKt = (4.*Kt_x.*ux.^(Chi_x+3).*(Chi_x+1).*(Bt_x.*Chi_x+Chi_x+2.*Bt_x+1))...
            ./(Fs_x.*(Bt_x+1).^2.*(Chi_x+2).*(Chi_x+3).*((Fs_x.*(Bt_x+1).*(Chi_x+2))...
            ./(Kt_x.*((Bt_x+1).*Chi_x+2.*Bt_x+1))).^Chi_x);
        
        dD1dChi = -(((4.*Kt_x.^2.*ux.^(Chi_x+3).*Bt_x+4.*Kt_x.^2.*ux.^(Chi_x+3)).*Chi_x.^3 ...
            +(24.*Kt_x.^2.*ux.^(Chi_x+3).*Bt_x+20.*Kt_x.^2.*ux.^(Chi_x+3)).*Chi_x.^2+(44.*Kt_x.^2.*ux.^(Chi_x+3).*Bt_x ...
            +28.*Kt_x.^2.*ux.^(Chi_x+3)).*Chi_x+24.*Kt_x.^2.*ux.^(Chi_x+3).*Bt_x+12.*Kt_x.^2.*ux.^(Chi_x+3)) ...
            .*log(((Fs_x.*Bt_x+Fs_x).*Chi_x+2.*Fs_x.*Bt_x+2.*Fs_x)./((Kt_x.*Bt_x+Kt_x).*Chi_x+2.*Kt_x.*Bt_x+Kt_x)) ...
            +(-4.*Kt_x.^2.*ux.^(Chi_x+3).*log(ux).*Bt_x-4.*Kt_x.^2.*ux.^(Chi_x+3).*log(ux)).*Chi_x.^3 ...
            +(ux.^Chi_x.*(4.*Kt_x.^2.*ux.^3-24.*Kt_x.^2.*ux.^3.*log(ux)).*Bt_x-20.*Kt_x.^2.*ux.^(Chi_x+3).*log(ux)).*Chi_x.^2 ...
            +(ux.^Chi_x.*(8.*Kt_x.^2.*ux.^3-44.*Kt_x.^2.*ux.^3.*log(ux)).*Bt_x ...
            +ux.^Chi_x.*(-28.*Kt_x.^2.*ux.^3.*log(ux)-8.*Kt_x.^2.*ux.^3)).*Chi_x ...
            +ux.^Chi_x.*(-24.*Kt_x.^2.*ux.^3.*log(ux)-4.*Kt_x.^2.*ux.^3).*Bt_x ...
            +ux.^Chi_x.*(-12.*Kt_x.^2.*ux.^3.*log(ux)-8.*Kt_x.^2.*ux.^3)) ...
            ./(((Fs_x.*Bt_x.^2+2.*Fs_x.*Bt_x+Fs_x).*Chi_x.^4+(10.*Fs_x.*Bt_x.^2+20.*Fs_x.*Bt_x+10.*Fs_x).*Chi_x.^3 ...
            +(37.*Fs_x.*Bt_x.^2+74.*Fs_x.*Bt_x+37.*Fs_x).*Chi_x.^2+(60.*Fs_x.*Bt_x.^2+120.*Fs_x.*Bt_x+60.*Fs_x).*Chi_x ...
            +36.*Fs_x.*Bt_x.^2+72.*Fs_x.*Bt_x+36.*Fs_x)...
            .*(((Fs_x.*Bt_x+Fs_x).*Chi_x+2.*Fs_x.*Bt_x+2.*Fs_x)./((Kt_x.*Bt_x+Kt_x).*Chi_x+2.*Kt_x.*Bt_x+Kt_x)).^Chi_x);
        
        dD1dBt = -(4.*Kt_x.^2.*ux.^(Chi_x+3).*Bt_x.*(Chi_x+1))./(Fs_x.*(Bt_x+1).^3.*(Chi_x+2).*(Chi_x+3)...
            .*((Fs_x.*(Bt_x+1).*(Chi_x+2))./(Kt_x.*((Bt_x+1).*Chi_x+2.*Bt_x+1))).^Chi_x);
    
        dD2dFs = (4.*(Kt_x.*ux.*Bt_x.^2.*Chi_x.^3-2.*Fs_x.*Bt_x.^2.*Chi_x.^3+2.*Kt_x.*ux.*Bt_x.*Chi_x.^3-4.*Fs_x.*Bt_x.*Chi_x.^3 ...
                +Kt_x.*ux.*Chi_x.^3-2.*Fs_x.*Chi_x.^3+7.*Kt_x.*ux.*Bt_x.^2.*Chi_x.^2-14.*Fs_x.*Bt_x.^2.*Chi_x.^2+12.*Kt_x.*ux.*Bt_x.*Chi_x.^2 ...
                -24.*Fs_x.*Bt_x.*Chi_x.^2+5.*Kt_x.*ux.*Chi_x.^2-10.*Fs_x.*Chi_x.^2+16.*Kt_x.*ux.*Bt_x.^2.*Chi_x-32.*Fs_x.*Bt_x.^2.*Chi_x ...
                +22.*Kt_x.*ux.*Bt_x.*Chi_x-48.*Fs_x.*Bt_x.*Chi_x+7.*Kt_x.*ux.*Chi_x-16.*Fs_x.*Chi_x+12.*Kt_x.*ux.*Bt_x.^2 ...
                -24.*Fs_x.*Bt_x.^2+12.*Kt_x.*ux.*Bt_x-32.*Fs_x.*Bt_x+3.*Kt_x.*ux-8.*Fs_x)) ...
                ./(Kt_x.*(Chi_x+3).*(Bt_x.*Chi_x+Chi_x+2.*Bt_x+1).^2);
        
        dD2dKt = (4.*Fs_x.^2.*(Bt_x+1).*(Chi_x+2).^2.*(Bt_x.*Chi_x+Chi_x+3.*Bt_x+1)) ...
                ./(Kt_x.^2.*(Chi_x+3).*(Bt_x.*Chi_x+Chi_x+2.*Bt_x+1).^2);
        
        dD2dChi = (8.*Fs_x.^2.*(Bt_x+1).*(Chi_x+2).*(2.*Bt_x.*Chi_x+Chi_x+5.*Bt_x+1)) ...
                    ./(Kt_x.*(Chi_x+3).^2.*(Bt_x.*Chi_x+Chi_x+2.*Bt_x+1).^3);
        
        dD2dBt = (8.*Fs_x.^2.*Bt_x.*(Chi_x+2).^2)...
                ./(Kt_x.*(Chi_x+3).*(Bt_x.*Chi_x+Chi_x+2.*Bt_x+1).^3);
        
        varargout{2}(:, :, 1) = [opt.x.T{1}(:, 1).*( dD1dFs.*(ux<PhiMx_x) + dD2dFs.*(ux>=PhiMx_x) ), ... %Fs
                                 opt.x.T{2}(:, 2).*( dD1dKt.*(ux<PhiMx_x) + dD2dKt.*(ux>=PhiMx_x) ), ... %Kt
                                 opt.x.T{3}(:, 3).*( dD1dChi.*(ux<PhiMx_x) + dD2dChi.*(ux>=PhiMx_x) ), ... %Chi
                                 opt.x.T{4}(:, 4).*( dD1dBt.*(ux<PhiMx_x) + dD2dBt.*(ux>=PhiMx_x) ), ... %Bt
                                 zeros(size(uxyn(:, 1)))]; %derivative w.r.t. Kn
        
        varargout{2}(:, :, 1) = varargout{2}(:, :, 1).*dlparsdpars';
                
        
        % Iwan 4 - Y
        
        dD1dFs = -(4.*Kt_y.^2.*uy.^(Chi_y+3).*(Chi_y+1).^2.*(Bt_y.*Chi_y+Chi_y+2.*Bt_y+1))...
            ./(Fs_y.^2.*(Bt_y+1).^2.*(Chi_y+2).^2.*(Chi_y+3).*((Fs_y.*(Bt_y+1).*(Chi_y+2))...
            ./(Kt_y.*((Bt_y+1).*Chi_y+2.*Bt_y+1))).^Chi_y);
        
        dD1dKt = (4.*Kt_y.*uy.^(Chi_y+3).*(Chi_y+1).*(Bt_y.*Chi_y+Chi_y+2.*Bt_y+1))...
            ./(Fs_y.*(Bt_y+1).^2.*(Chi_y+2).*(Chi_y+3).*((Fs_y.*(Bt_y+1).*(Chi_y+2))...
            ./(Kt_y.*((Bt_y+1).*Chi_y+2.*Bt_y+1))).^Chi_y);
        
        dD1dChi = -(((4.*Kt_y.^2.*uy.^(Chi_y+3).*Bt_y+4.*Kt_y.^2.*uy.^(Chi_y+3)).*Chi_y.^3 ...
            +(24.*Kt_y.^2.*uy.^(Chi_y+3).*Bt_y+20.*Kt_y.^2.*uy.^(Chi_y+3)).*Chi_y.^2+(44.*Kt_y.^2.*uy.^(Chi_y+3).*Bt_y ...
            +28.*Kt_y.^2.*uy.^(Chi_y+3)).*Chi_y+24.*Kt_y.^2.*uy.^(Chi_y+3).*Bt_y+12.*Kt_y.^2.*uy.^(Chi_y+3)) ...
            .*log(((Fs_y.*Bt_y+Fs_y).*Chi_y+2.*Fs_y.*Bt_y+2.*Fs_y)./((Kt_y.*Bt_y+Kt_y).*Chi_y+2.*Kt_y.*Bt_y+Kt_y)) ...
            +(-4.*Kt_y.^2.*uy.^(Chi_y+3).*log(uy).*Bt_y-4.*Kt_y.^2.*uy.^(Chi_y+3).*log(uy)).*Chi_y.^3 ...
            +(uy.^Chi_y.*(4.*Kt_y.^2.*uy.^3-24.*Kt_y.^2.*uy.^3.*log(uy)).*Bt_y-20.*Kt_y.^2.*uy.^(Chi_y+3).*log(uy)).*Chi_y.^2 ...
            +(uy.^Chi_y.*(8.*Kt_y.^2.*uy.^3-44.*Kt_y.^2.*uy.^3.*log(uy)).*Bt_y ...
            +uy.^Chi_y.*(-28.*Kt_y.^2.*uy.^3.*log(uy)-8.*Kt_y.^2.*uy.^3)).*Chi_y ...
            +uy.^Chi_y.*(-24.*Kt_y.^2.*uy.^3.*log(uy)-4.*Kt_y.^2.*uy.^3).*Bt_y ...
            +uy.^Chi_y.*(-12.*Kt_y.^2.*uy.^3.*log(uy)-8.*Kt_y.^2.*uy.^3)) ...
            ./(((Fs_y.*Bt_y.^2+2.*Fs_y.*Bt_y+Fs_y).*Chi_y.^4+(10.*Fs_y.*Bt_y.^2+20.*Fs_y.*Bt_y+10.*Fs_y).*Chi_y.^3 ...
            +(37.*Fs_y.*Bt_y.^2+74.*Fs_y.*Bt_y+37.*Fs_y).*Chi_y.^2+(60.*Fs_y.*Bt_y.^2+120.*Fs_y.*Bt_y+60.*Fs_y).*Chi_y ...
            +36.*Fs_y.*Bt_y.^2+72.*Fs_y.*Bt_y+36.*Fs_y)...
            .*(((Fs_y.*Bt_y+Fs_y).*Chi_y+2.*Fs_y.*Bt_y+2.*Fs_y)./((Kt_y.*Bt_y+Kt_y).*Chi_y+2.*Kt_y.*Bt_y+Kt_y)).^Chi_y);
        
        dD1dBt = -(4.*Kt_y.^2.*uy.^(Chi_y+3).*Bt_y.*(Chi_y+1))./(Fs_y.*(Bt_y+1).^3.*(Chi_y+2).*(Chi_y+3)...
            .*((Fs_y.*(Bt_y+1).*(Chi_y+2))./(Kt_y.*((Bt_y+1).*Chi_y+2.*Bt_y+1))).^Chi_y);
    
        dD2dFs = (4.*(Kt_y.*uy.*Bt_y.^2.*Chi_y.^3-2.*Fs_y.*Bt_y.^2.*Chi_y.^3+2.*Kt_y.*uy.*Bt_y.*Chi_y.^3-4.*Fs_y.*Bt_y.*Chi_y.^3 ...
                +Kt_y.*uy.*Chi_y.^3-2.*Fs_y.*Chi_y.^3+7.*Kt_y.*uy.*Bt_y.^2.*Chi_y.^2-14.*Fs_y.*Bt_y.^2.*Chi_y.^2+12.*Kt_y.*uy.*Bt_y.*Chi_y.^2 ...
                -24.*Fs_y.*Bt_y.*Chi_y.^2+5.*Kt_y.*uy.*Chi_y.^2-10.*Fs_y.*Chi_y.^2+16.*Kt_y.*uy.*Bt_y.^2.*Chi_y-32.*Fs_y.*Bt_y.^2.*Chi_y ...
                +22.*Kt_y.*uy.*Bt_y.*Chi_y-48.*Fs_y.*Bt_y.*Chi_y+7.*Kt_y.*uy.*Chi_y-16.*Fs_y.*Chi_y+12.*Kt_y.*uy.*Bt_y.^2 ...
                -24.*Fs_y.*Bt_y.^2+12.*Kt_y.*uy.*Bt_y-32.*Fs_y.*Bt_y+3.*Kt_y.*uy-8.*Fs_y)) ...
                ./(Kt_y.*(Chi_y+3).*(Bt_y.*Chi_y+Chi_y+2.*Bt_y+1).^2);
        
        dD2dKt = (4.*Fs_y.^2.*(Bt_y+1).*(Chi_y+2).^2.*(Bt_y.*Chi_y+Chi_y+3.*Bt_y+1)) ...
                ./(Kt_y.^2.*(Chi_y+3).*(Bt_y.*Chi_y+Chi_y+2.*Bt_y+1).^2);
        
        dD2dChi = (8.*Fs_y.^2.*(Bt_y+1).*(Chi_y+2).*(2.*Bt_y.*Chi_y+Chi_y+5.*Bt_y+1)) ...
                    ./(Kt_y.*(Chi_y+3).^2.*(Bt_y.*Chi_y+Chi_y+2.*Bt_y+1).^3);
        
        dD2dBt = (8.*Fs_y.^2.*Bt_y.*(Chi_y+2).^2)...
                ./(Kt_y.*(Chi_y+3).*(Bt_y.*Chi_y+Chi_y+2.*Bt_y+1).^3);
        
        varargout{2}(:, :, 2) = [opt.y.T{1}(:, 1).*( dD1dFs.*(uy<PhiMx_y) + dD2dFs.*(uy>=PhiMx_y) ), ... %Fs
                                 opt.y.T{2}(:, 2).*( dD1dKt.*(uy<PhiMx_y) + dD2dKt.*(uy>=PhiMx_y) ), ... %Kt
                                 opt.y.T{3}(:, 3).*( dD1dChi.*(uy<PhiMx_y) + dD2dChi.*(uy>=PhiMx_y) ), ... %Chi
                                 opt.y.T{4}(:, 4).*( dD1dBt.*(uy<PhiMx_y) + dD2dBt.*(uy>=PhiMx_y) ), ... %Bt
                                 zeros(size(uxyn(:, 1)))]; %derivative w.r.t. Kn
        
        varargout{2}(:, :, 2) = varargout{2}(:, :, 2).*dlparsdpars';
        
        % Spring - N: No dissipation
        varargout{2}(:, :, 3) = zeros(size(Dxyn,1), length(lpars));
    end
end
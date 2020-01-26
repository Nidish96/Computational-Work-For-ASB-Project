function [Fxyn, varargout] = SPRING_IWAN4(uxyn, pars, opt)

    lpars = pars(:);
    lpars(opt.lspci) = 10.^(lpars(opt.lspci));
    dlparsdpars = ones(size(lpars));
    dlparsdpars(opt.lspci) = lpars(opt.lspci).*log(10);    
    
    %% Forces    
    Fxyn = zeros(size(uxyn));
    
    % Iwan 4 - X
    Fs = opt.x.T{1}*lpars;
    Kt = opt.x.T{2}*lpars;
    Chi = opt.x.T{3}*lpars;
    Bt = opt.x.T{4}*lpars;
    
    PhiMx = Fs.*(1+Bt)./(Kt.*(Bt+(Chi+1)./(Chi+2)));
    CoefE = (Kt.*(Kt.*(Bt+(Chi+1)./(Chi+2))./(Fs.*(1+Bt))).^(1+Chi))./((1+Bt).*(Chi+2));    
    
    sux = sign(real(uxyn(:,1)));
    ux = sux.*uxyn(:,1);
	Fxyn(:,1) = ((Kt.*ux - CoefE.*ux.^(2+Chi)).*(ux<PhiMx) + (Fs).*(ux>=PhiMx)).*sux;
    
    % Iwan 4 - Y
    fs = opt.y.T{1}*lpars;
    kt = opt.y.T{2}*lpars;
    chi = opt.y.T{3}*lpars;
    bt = opt.y.T{4}*lpars;
    
    phiMx = fs.*(1+bt)./(kt.*(bt+(chi+1)./(chi+2)));
    coefE = (kt.*(kt.*(bt+(chi+1)./(chi+2))./(fs.*(1+bt))).^(1+chi))./((1+bt).*(chi+2));
   
    suy = sign(real(uxyn(:,2)));
    uy = suy.*uxyn(:,2);
	Fxyn(:,2) = ((kt.*uy - coefE.*uy.^(2+chi)).*(uy<phiMx) + (fs).*(uy>=phiMx)).*suy;
    
    % Linear Spring - N
    Kn = opt.n.T{1}*lpars;
    Fxyn(:, 3) = Kn.*uxyn(:,3);
    %% Jacobians
    if nargout>=2
        varargout{1} = [(Kt-CoefE.*(2+Chi).*ux.^(1+Chi)).*(ux<PhiMx), (kt-coefE.*(2+chi).*uy.^(1+chi)).*(uy<phiMx), Kn];
    end
    if nargout>=3
        varargout{2} = zeros(size(Fxyn,1), length(lpars), 3);
        
        % Iwan 4 - X
        t1 = log(Kt.*(Bt+(Chi+1)./(Chi+2))./(Fs.*(Bt+1)));
        t2 = Bt./((Bt+1).*Chi+2*Bt+1);
        dCdp = CoefE.*[-(1+Chi)./Fs, ...
                        (2+Chi)./Kt, ...
                        t1-t2, ...
                        (-Bt./((Bt+1).*(Bt+(1+Chi)./(2+Chi))))];            
        tmp = [-dCdp(:,1).*ux.^(2+Chi), ...
                ux-dCdp(:,2).*ux.^(2+Chi), ...
              -(dCdp(:,3)+CoefE.*log(ux)).*ux.^(Chi+2), ...
               -dCdp(:,4).*ux.^(2+Chi)].*(ux<PhiMx) + ...
             [1 0 0 0].*(ux>=PhiMx);
        tmp = sux.*tmp;
        for i=1:size(tmp,2)
            varargout{2}(:, :, 1) = varargout{2}(:, :, 1) + diag(tmp(:,i))*opt.x.T{i};
        end
        varargout{2}(:,:,1) = varargout{2}(:,:,1).*dlparsdpars';  % Check
        
        % Iwan 4 - Y
        t1 = log(kt.*(bt+(chi+1)./(chi+2))./(fs.*(bt+1)));
        t2 = bt./((bt+1).*chi+2*bt+1);
        dCdp = coefE.*[-(1+chi)./fs, ...
                        (2+chi)./kt, ...
                        t1-t2, ...
                        -bt./((bt+1).*(bt+(1+chi)./(2+chi)))];
            
        tmp = [-dCdp(:,1).*uy.^(2+chi), ...
                uy-dCdp(:,2).*uy.^(2+chi), ...
              -(dCdp(:,3)+coefE.*log(uy)).*uy.^(chi+2), ...
               -dCdp(:,4).*uy.^(2+chi)].*(uy<phiMx) + ...
             [1 0 0 0].*(uy>=phiMx);
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
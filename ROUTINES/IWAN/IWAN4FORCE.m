function [F, varargout] = IWAN4FORCE(ut, pars, totcol, occol)
    Fs = pars(:, occol(1));  
    Kt = pars(:, occol(2));  
    Chi = pars(:, occol(3)); 
    Bt = pars(:, occol(4));
    
    PhiMx = Fs.*(1+Bt)./(Kt.*(Bt+(Chi+1)./(Chi+2)));
    CoefE = (Kt.*(Kt.*(Bt+(Chi+1)./(Chi+2))./(Fs.*(1+Bt))).^(1+ ...
                                                      Chi))./((1+Bt).*(Chi+2));    
    
    sut = sign(real(ut));
    ut = sut.*ut;
    F = ((Kt.*ut - CoefE.*ut.^(2+Chi)).*(ut<PhiMx) + (Fs).*(ut>=PhiMx)).*sut;
    
    if nargout>=2  % Derivatives w.r.t displacements required
        varargout{1} = (Kt - CoefE.*(2+Chi).*ut.^(1+Chi)).*(ut<PhiMx);
        if nargout>=3  % Derivatives w.r.t parameters required
            t1 = log(Kt.*(Bt+(Chi+1)./(Chi+2))./(Fs.*(Bt+1)));
            t2 = Bt./((Bt+1).*Chi+2*Bt+1);
            dCdp = CoefE.*[-(1+Chi)./Fs, ...
                           (2+Chi)./Kt, ...
                           t1-t2, ...
                           (-Bt./((Bt+1).*(Bt+(1+Chi)./(2+Chi))))];
            
            tmp = [-dCdp(:,1).*ut.^(2+Chi), ...
                            ut-dCdp(:,2).*ut.^(2+Chi), ...
                            -(dCdp(:,3)+CoefE.*log(ut)).*ut.^(Chi+2), ...
                            -dCdp(:,4).*ut.^(2+Chi)].*(ut<PhiMx) + ...
                      [1 0 0 0].*(ut>=PhiMx);
            tmp = sut.*tmp;
            
            varargout{2} = zeros(size(tmp,1),totcol);
            varargout{2}(:, occol) = tmp;
        end
    end
end
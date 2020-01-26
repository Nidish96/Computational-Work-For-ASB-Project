function [F, dFdu] = JENK2FORCE(ut, pars)
    Fs = pars(:, 1);
    Kt = pars(:, 2);
    PhiMx = Fs./Kt;
    
    F = (Kt.*ut).*(abs(ut)<PhiMx) + ...
         (Fs.*sign(ut)).*(abs(ut)>=PhiMx);
     i = length(find(abs(ut)>PhiMx));
     disp([num2str(i) ' Slipped.'])
    
    dFdu = (Kt).*(abs(ut)<PhiMx);
    
    F = real(F);
    dFdu = real(dFdu);
end
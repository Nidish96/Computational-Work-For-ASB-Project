function [D] = IWAN4DISS_M(ut, pars)
    R = pars(:, 1);
    S = pars(:, 2);
    Chi = pars(:, 3);
    PhiMx = pars(:, 4);
    
    sut = sign(ut);
    ut = abs(ut);
    
    D = (4*R.*ut.^(Chi+3)./((Chi+2).*(Chi+3))).*(ut<PhiMx) + ...
        (4*R.*PhiMx.^(Chi+2).*(ut./(Chi+2)-PhiMx./(Chi+3))+ ...
        4*S.*PhiMx.*(ut-PhiMx)).*(ut>=PhiMx);
end
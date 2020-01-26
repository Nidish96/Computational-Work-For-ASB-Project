function [Werr, dWerrdp] = LOWAMPQSMA_WERRFUN(RESFUN, pars, alpha, Wexp, X0, Xstat, VmO, opt)    
try
    Xstat = fsolve(@(X) RESFUN([X; 0], pars, VmO'), Xstat, opt);
    [~,dRdX,~,dRdp] = RESFUN([Xstat; 0], pars, VmO');
    dXsdp = -dRdX\dRdp;

    X = fsolve(@(X) RESFUN([X; alpha], pars, VmO'), X0, opt);
    [~,dRdX,~,dRdp] = RESFUN([X; alpha], pars, VmO');
    dXdp = -dRdX\dRdp;
catch me
    disp('peekaboo');
end
    
    q = abs(VmO*(X-Xstat));
    W = sqrt(alpha/q);
    
    dWdp = -alpha/(2*W*q^2)*VmO*(dXdp-dXsdp);
    
    Werr = (W-Wexp)^2;
    dWerrdp = 2*(W-Wexp)*dWdp;
end
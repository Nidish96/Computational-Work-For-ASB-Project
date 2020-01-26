function [ErrsWZD, dErrsdp, BB] = WJMODEL_BBFUN(pars, mdi, expdat, K, M, X0, Fv, L, Txyn, Qxyn, CFUN, Npatches, Nqp, opt)
%WJMODEL_BBFUN returns the backbone

    % Prestress
    [Xstat, eflg] = fsolve(@(X) WJRES_STR([X; 0.0], pars, K, Fv, Fv, L, Txyn, Qxyn, CFUN, Npatches), X0, opt);
    if eflg<0
        disp(pars)
        error('Non-convergent!');
    end
    [~,dRstat,~,dRdpstat] = WJRES_STR([Xstat; 0.0], pars, K, Fv, Fv, L, Txyn, Qxyn, CFUN, Npatches);
    dXdpstat = -dRstat\dRdpstat;
	prev.uxyn = [Qxyn(1:3:end,:)*Xstat, Qxyn(2:3:end,:)*Xstat, Qxyn(3:3:end,:)*Xstat*0];
    [prev.Fxyn, dFxynstat, prev.dFxyndpar] = CFUN(prev.uxyn, pars);
	if size(dFxynstat,3)==1
       tmp = dFxynstat;
       dFxynstat = zeros(size(dFxynstat,1),size(dFxynstat,2),3);
       dFxynstat(:,1,1) = tmp(:,1);
       dFxynstat(:,2,2) = tmp(:,2);
       dFxynstat(:,3,3) = tmp(:,3); % Np x uxyn x fxyn
    end
	dFxynstat = num2cell(permute(dFxynstat,[3,2,1]), [2,1]);
    prev.dFxyndXdpar = blkdiag(dFxynstat{:})*Qxyn*dXdpstat;
    
%     prev.dFxyndXdpar = reshape(dFxynstat', Npatches*3, 1).*Qxyn*dXdpstat;
    
%     prev.uxyn = prev.uxyn*0;
%     prev.Fxyn = prev.Fxyn*0;
%     prev.dFxyndpar = prev.dFxyndpar*0;
%     prev.dFxyndXdpar = prev.dFxyndXdpar*0;
    %% Modal Analysis
    [Vst, Wst] = eigs(dRstat, M, 10, 'SM');
    Wst = sqrt(diag(Wst));
    Vst = Vst./sqrt(diag(Vst'*M*Vst)');

    Vm = Vst(:, mdi);
    Wm = Wst(mdi);
    %% RQNM
    BB.Q = expdat.Q;
    BB.W = zeros(size(BB.Q));
    BB.Z = zeros(size(BB.Q));
    BB.D = zeros(size(BB.Q));
    
    dWdp = zeros(1, length(pars));
    dDdp = zeros(1, length(pars));
    dZdp = zeros(1, length(pars));
    
    ErrsWZD = zeros(3, 1);
    dErrsdp = zeros(3, length(pars));
    
    % LGL Quadrature points
    [xi, wi] = LGLWT(Nqp, 0, 1);
    qr = zeros(size(xi));
    
    %alternative trapezoids: - this showed very slow convergence and may
    %requires > 400 points to converge. 
%     xi = linspace(0, 1, Nqp+1)';
%     wi = 1/Nqp/2*[1; 2*ones(Nqp-1, 1);1];
%     qr = zeros(size(xi));
    
    % Initialization for accumulators
    Q0 = 0;
    Di = 0;
    dDdpi = zeros(1, length(pars));
    Vmi = Vm;
    Wmi = Wm;
    
    rmp = 0;  % Number of points to remove
    for i=1:length(BB.Q)
        if i>1
            Q0 = BB.Q(i-1);
        end
        
        %Set things to zero for restarting all quadature integartion from
        %0:
%         BB.D(i) = 0;
%         dDdp = zeros(1, length(pars));
%         Q0 = 0;
%         Di = 0;
        
        for k=2:(Nqp+1)  % Evaluate Nqp points for quadrature
            qr(k) = Q0 + (BB.Q(i)-Q0)*xi(k);
            X0 = [Xstat+Vmi*qr(k); Wmi^2];
            Xd = fsolve(@(X) WJRES_RQNM([X; qr(k)], pars, prev, K, M, Xstat, dXdpstat, ...
                Fv, L, Txyn, Qxyn, CFUN, Npatches), X0, opt);
            [~, dRdX, ~, dRdp] = WJRES_RQNM([Xd; qr(k)], pars, prev, K, M, Xstat, dXdpstat, Fv, ...
                L, Txyn, Qxyn, CFUN, Npatches);
            dXdp = -dRdX\dRdp;
            
            Di = Di + (BB.Q(i)-Q0)*wi(k)*qr(k)*Xd(end);
            dDdpi = dDdpi + (BB.Q(i)-Q0)*wi(k)*qr(k)*dXdp(end,:);
        end
        Vmi = (Xd(1:end-1)-Xstat)/BB.Q(i);
        Wmi = sqrt(Xd(end));
        
        BB.D(i) = -4*BB.Q(i)^2*Xd(end) + 8*Di;
        dDdp = -4*BB.Q(i)^2*dXdp(end,:) + 8*dDdpi;
        
        BB.W(i) = sqrt(Xd(end));
        BB.Z(i) = BB.D(i)/(2*pi*(BB.Q(i)*BB.W(i))^2) + expdat.Z(1);
        
        dWdp = dXdp(end,:)/(2*BB.W(i));
        dZdp = (dDdp - 2*BB.D(i)/BB.W(i)*dWdp)/(2*pi*(BB.Q(i)*BB.W(i))^2);
        
        if isreal(BB.W(i))
            ErrsWZD(1) = ErrsWZD(1) + ((BB.W(i)-expdat.W(i))/expdat.W(i))^2;
            dErrsdp(1,:) = dErrsdp(1,:) + 2*((BB.W(i)-expdat.W(i))/expdat.W(i)^2)*dWdp;
        else
            ErrsWZD(1) = ErrsWZD(1) + 1;
            dErrsdp(1,:) = dErrsdp(1,:) + 0;
        end
        
%         ErrsWZD(2) = ErrsWZD(2) + ((BB.Z(i)-expdat.Z(i))/expdat.Z(i))^2;
%         ErrsWZD(3) = ErrsWZD(3) + ((BB.D(i)-expdat.D(i))/expdat.D(i))^2;
%         dErrsdp(2,:) = dErrsdp(2,:) + 2*((BB.Z(i)-expdat.Z(i))/expdat.Z(i)^2)*dZdp;
%         dErrsdp(3,:) = dErrsdp(3,:) + 2*((BB.D(i)-expdat.D(i))/expdat.D(i)^2)*dDdp;
        
        if BB.Z(i)>0
            ErrsWZD(2) = ErrsWZD(2) + ((log(BB.Z(i))-log(expdat.Z(i)))/log(expdat.Z(i)))^2;
            ErrsWZD(3) = ErrsWZD(3) + ((log(BB.D(i))-log(expdat.D(i)))/log(expdat.D(i)))^2;
            dErrsdp(2,:) = dErrsdp(2,:) + 2/BB.Z(i)*((log(BB.Z(i))-log(expdat.Z(i)))/log(expdat.Z(i))^2)*dZdp;
            dErrsdp(3,:) = dErrsdp(3,:) + 2/BB.Z(i)*((log(BB.D(i))-log(expdat.D(i)))/log(expdat.D(i))^2)*dDdp;
        else
            ErrsWZD(2) = ErrsWZD(2) + 1;
            ErrsWZD(3) = ErrsWZD(3) + 1;
            dErrsdp(2,:) = dErrsdp(2,:) + 0;
            dErrsdp(3,:) = dErrsdp(3,:) + 0;
        end
%         fprintf('%d ',i);
    end
%     fprintf('\n');
    ErrsWZD = ErrsWZD/length(BB.Q);
    dErrsdp = dErrsdp/length(BB.Q);
    
    ErrsWZD = ErrsWZD([1 2],:);
    dErrsdp = dErrsdp([1 2],:);

    % Log-Scale Errors
    dErrsdp = (dErrsdp./ErrsWZD)/log(10);
    ErrsWZD = log10(ErrsWZD);
    
    if ~isreal(ErrsWZD)
        disp('complex');
    end
end

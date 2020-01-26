function [R, dRdX, dRdq, varargout] = WJRES_RQNM(X, pars, prev, K, M, Xstat, dXsdp, Fs, L, Txyn, Qxyn, Fxynfun, Nps, varargin)
%WJRESFUN Residual Function for Whole Joint Models
% USAGE:
%  [R, dRdX, dRda] = WJRESFUN(X, K, Fd, Fs, L, Fxfun, Fyfun, Fnfun);
% INPUTS:
%   X       :
%   K       :
% OUTPUTS:
%   R
%   dRdX
%   dRda
    if nargin>=14
        RCvec = varargin{1};
        RX = RCvec.*X;
    else
        RCvec = ones(size(X));
        RX = X;
    end
    
    q = RX(end);
    lam = RX(end-1);
    
    uxyn = [Qxyn(1:3:end,:)*RX(1:end-2) Qxyn(2:3:end,:)*RX(1:end-2) Qxyn(3:3:end,:)*RX(1:end-2)];
    if nargout<4
        [Fxyn, dFxyn] = Fxynfun(uxyn-prev.uxyn, pars);
        Fxyn = Fxyn + prev.Fxyn;
    else
        [Fxyn, dFxyn, dFxyndpar] = Fxynfun(uxyn-prev.uxyn, pars);
        Fxyn = Fxyn + prev.Fxyn;
        dFxyndpar = dFxyndpar + prev.dFxyndpar;
    end
    if size(dFxyn,3)==1
       tmp = dFxyn;
       dFxyn = zeros(size(dFxyn,1),size(dFxyn,2),3);
       dFxyn(:,1,1) = tmp(:,1);
       dFxyn(:,2,2) = tmp(:,2);
       dFxyn(:,3,3) = tmp(:,3); % Np x uxyn x fxyn
    end
	dFxyn = num2cell(permute(dFxyn,[2,3,1]), [2,1]);
    
    Ng = size(L,1) - Nps*6;
    Fnl = Txyn*reshape(Fxyn', Nps*3, 1);
%     dFnl = Txyn*diag(reshape(dFxyn', Nps*3, 1))*Qxyn;
    dFnl = Txyn*blkdiag(dFxyn{:})*Qxyn;
    
    R = [K*RX(1:end-2) + Fnl - Fs - lam*M*(RX(1:end-2)-Xstat);
        (RX(1:end-2)-Xstat)'*M*(RX(1:end-2)-Xstat)-q^2];
    dRdX = [(K - lam*M + dFnl).*RCvec(1:end-2), -M*(RX(1:end-2)-Xstat);
        2*(RX(1:end-2)-Xstat)'*M, 0];
    dRdq = [zeros(length(RX)-2,1); -2*q];
    
	if nargin>=15
        R = varargin{2}.*R;
        dRdX = varargin{2}.*dRdX;
        dRdq = varargin{2}.*dRdq;
    end
    if nargout>=4
        varargout{1} = zeros(Nps*3, length(pars));
        for i=1:length(pars)
            varargout{1}(:, i) = reshape(squeeze(dFxyndpar(:,i,:))', Nps*3, 1);
        end
%         varargout{1} = reshape(permute(dFxyndpar, [3 1 2]), 3*Nps, length(pars))-(reshape(dFxyn', Nps*3, 1)).*Qxyn*dXsdp+prev.dFxyndXdpar;
        varargout{1} = reshape(permute(dFxyndpar, [3 1 2]), 3*Nps, length(pars))-blkdiag(dFxyn{:})*Qxyn*dXsdp+prev.dFxyndXdpar;
        varargout{1} = [Txyn*varargout{1}+lam*M*dXsdp; -2*(RX(1:end-2)-Xstat)'*M*dXsdp];
    end
end
function [R, dRdX, dRda, varargout] = WJRES_STR(X, pars, K, Fd, Fs, L, Txyn, Qxyn, Fxynfun, Nps, varargin)
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
    if nargin>=13
        RCvec = varargin{1};
        RX = RCvec.*X;
    else
        RCvec = ones(size(X));
        RX = X;
    end
    
    alpha = RX(end);
    uxyn = [Qxyn(1:3:end,:)*RX(1:end-1) Qxyn(2:3:end,:)*RX(1:end-1) Qxyn(3:3:end,:)*RX(1:end-1)];
    if nargout<4
        [Fxyn, dFxyn] = Fxynfun(uxyn, pars);
    else        
        [Fxyn, dFxyn, dFxyndpar] = Fxynfun(uxyn, pars);
    end
	if size(dFxyn,3)==1
       tmp = dFxyn;
       dFxyn = zeros(size(dFxyn,1),size(dFxyn,2),3);
       dFxyn(:,1,1) = tmp(:,1);
       dFxyn(:,2,2) = tmp(:,2);
       dFxyn(:,3,3) = tmp(:,3); % Np x uxyn x fxyn
    end
	dFxyn = num2cell(permute(dFxyn,[3,2,1]), [2,1]);
    
    Ng = size(L,1) - Nps*6;
    Fnl = Txyn*reshape(Fxyn', Nps*3, 1);
%     dFnl = Txyn*diag(reshape(dFxyn', Nps*3, 1))*Qxyn;
	dFnl = Txyn*blkdiag(dFxyn{:})*Qxyn;
    
    R = K*RX(1:end-1) + Fnl - Fs - alpha*Fd;
    dRdX = K.*RCvec(1:end-1) + dFnl.*RCvec(1:end-1);
    dRda = -Fd*RCvec(end);
    
	if nargin>=14
        R = varargin{2}.*R;
        dRdX = varargin{2}.*dRdX;
        dRda = varargin{2}.*dRda;
    end
    if nargout>=4
        varargout{1} = Txyn*reshape(permute(dFxyndpar, [3, 1, 2]), Nps*3, length(pars));
    end
end

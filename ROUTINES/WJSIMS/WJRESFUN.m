function [R, dRdX, dRda, varargout] = WJRESFUN(X, pars, K, Fd, Fs, L, Txyn, Qxyn, Fxfun, Fyfun, Fnfun, Nps, varargin)
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
    if nargout<4
        [Fx, dFx] = Fxfun(Qxyn(1:3:end,:)*RX(1:end-1), pars);
        [Fy, dFy] = Fyfun(Qxyn(2:3:end,:)*RX(1:end-1), pars);
        [Fn, dFn] = Fnfun(Qxyn(3:3:end,:)*RX(1:end-1), pars);
    else
        [Fx, dFx, dFxdp] = Fxfun(Qxyn(1:3:end,:)*RX(1:end-1), pars);
        [Fy, dFy, dFydp] = Fyfun(Qxyn(2:3:end,:)*RX(1:end-1), pars);
        [Fn, dFn, dFndp] = Fnfun(Qxyn(3:3:end,:)*RX(1:end-1), pars);
    end
    
    Ng = size(L,1) - Nps*6;
    Fnl = Txyn*reshape([Fx Fy Fn]', Nps*3, 1);
    dFnl = Txyn*diag(reshape([dFx dFy dFn]', Nps*3, 1))*Qxyn;
    
    R = K*RX(1:end-1) + Fnl - Fs - alpha*Fd;
    dRdX = K.*RCvec(1:end-1) + dFnl.*RCvec(1:end-1);
    dRda = -Fd*RCvec(end);
    
	if nargin>=14
        R = varargin{2}.*R;
        dRdX = varargin{2}.*dRdX;
        dRda = varargin{2}.*dRda;
    end
    if nargout>=4
        varargout{1} = zeros(Nps*3, length(pars));
        for i=1:5
            varargout{1}(:, i) = reshape([dFxdp(:,i) dFydp(:,i) dFndp(:,i)]', Nps*3, 1);
        end
        varargout{1} = dRdX\(Txyn*varargout{1});
    end
end


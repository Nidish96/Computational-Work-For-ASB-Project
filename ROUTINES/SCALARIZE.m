function [varargout] = SCALARIZE(fun, X, varargin)
%SCALARIZE2 scalarizes the multi-objective optimization problem as a single
%objective problem.
%  USAGE:
%       [O, dOdp, dOdt] = SCALARIZE(fun, X, varargin)
%  INPUTS:
%       fun
%       X
%       opts (optional)
%  OUTPUTS:
%       O
%       dOdX
%       dOdp
    
    % Default Options
    opts = struct('method', 'weighted', 'gradients', true, 'nobj', 2, 'npar', 1, 'nvar', 1);
    if nargin==3
        nflds = fieldnames(varargin{1});
        for i=1:length(nflds)
            opts.(nflds{i}) = varargin{1}.(nflds{i});
        end
    end    
    % Evaluate all objectives
    if opts.gradients
        [Os, dOsdX] = fun(X(1:opts.nvar));
    else
        Os = fun(X(1:end-1));
        dOsdX = zeros(length(Os), length(X(1:opts.nvar)));
    end
    
    % Conduct scalarization
    pars = X((end-(opts.npar-1)):end);
    switch opts.method
        case 'weighted'
            if opts.npar~=opts.nobj-1
                error('Unrecognized parameters');
            end
            varargout{1} = [pars 1-sum(pars)]*Os(1:opts.nobj);
            varargout{2} = [pars 1-sum(pars)]*dOsdX(1:opts.nobj,:);
            varargout{3} = Os(1:opts.nobj)'*[eye(opts.npar); -ones(1,opts.npar)];
        case 'sphericalwgts'
            if opts.npar~=opts.nobj-1
                error('Unrecognized parameters');
            end            
            if ~isfield(opts, 'rpt')
                error('Source point not given for spherical boundary intersection');
            end
            rhat = ones(opts.nobj,1);
            drhatdp = zeros(opts.nobj, opts.npar);
            for i=1:opts.npar
                drhatdp(i, :) = drhatdp(i, :)*cos(pars(i)) - rhat(i)*[zeros(1,opts.npar-i-1), sin(pars(i)), zeros(1,opts.npar-i)];
                drhatdp((i+1):end, :) = drhatdp((i+1):end, :).*sin(pars(i)) + rhat((i+1):end).*[zeros(1,opts.npar-i-1), cos(pars(i)), zeros(1,opts.npar-i)];
                
                rhat(i) = rhat(i)*cos(pars(i));
                rhat((i+1):end) = rhat((i+1):end).*sin(pars(i));
            end
            drhatdp = (sum(rhat)*drhatdp-rhat.*sum(drhatdp,1))/sum(rhat)^2;
            rhat = rhat/sum(rhat);
            
            varargout{1} = rhat'*(Os(1:opts.nobj)-opts.rpt);
            varargout{2} = rhat'*dOsdX(1:opts.nobj,:);
            varargout{3} = (Os(1:opts.nobj)-opts.rpt)'*drhatdp;
        case 'sphericalbi'
            if opts.npar~=opts.nobj
                error('Unrecognized parameters');
            end
            if ~isfield(opts, 'rpt')
                error('Source point not given for spherical boundary intersection');
            end
            rhat = ones(opts.nobj,1);
            drhatdp = zeros(opts.nobj, opts.npar-1);
            for i=1:opts.npar-1  % Last Parameter is radial distance
                drhatdp(i, :) = drhatdp(i, :)*cos(pars(i+1)) - rhat(i)*[zeros(1,opts.npar-i-2), sin(pars(i+1)), zeros(1,opts.npar-i-1)];
                drhatdp((i+1):end, :) = drhatdp((i+1):end, :).*sin(pars(i+1)) + rhat((i+1):end).*[zeros(1,opts.npar-i-2), cos(pars(i+1)), zeros(1,opts.npar-i-1)];
                
                rhat(i) = rhat(i)*cos(pars(i+1));
                rhat((i+1):end) = rhat((i+1):end).*sin(pars(i+1));
            end
            drhatdp = (sum(rhat)*drhatdp-rhat.*sum(drhatdp,1))/sum(rhat)^2;
            rhat = rhat/sum(rhat);
            
%             O = opts.rpt + pars(end)*rhat - Os;
%             dOdX = -dOsdX;
%             dOdp = [pars(end)*drhatdp rhat];
            
            varargout{1} = [];
            varargout{2} = opts.rpt + pars(end-1)*rhat - Os;
            varargout{3} = [];
            varargout{4} = [-dOsdX rhat]';
            varargout{5} = [pars(end-1)*drhatdp rhat];
    end
end







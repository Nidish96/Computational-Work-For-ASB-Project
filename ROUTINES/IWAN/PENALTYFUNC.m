function [F,varargout] = PENALTYFUNC(un, pars, totcol, occol)
    Kn = pars(:,occol(1));
    F = max(Kn.*un, 0);
    
    if nargout>=2
        varargout{1} = Kn.*(F~=0);
        
        if nargout>=3
            tmp = zeros(length(un), totcol);
            tmp(:, occol) = un.*(F~=0);
            varargout{2} = tmp;
        end
    end
end
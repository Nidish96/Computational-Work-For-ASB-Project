function [Felem] = TRI2D_ELINT(func,V,no,nforces,nderivs)
%TRI2D_ELINT Integrate given scalar-valued function over whole
%domain and report nodal components
% USAGE:
%	[Vnodal] = TRI2D_ELINT(func,V,no);
% INPUTS:
%   func	: (1x1)-valued vectorized function (in global CS)
%       fval = func([x y]);
%       function could have multiple output scalars - see nouts
%   V 		: 3x2 vertex coordinates (global CS)
%   no		: (int) number of points for GQL - forced to be a
%   		  perfect square
%   nouts   : Number of output arguments of function
% OUTPUTS:
%   Fnodal	: noutsx1 cell of 3x1 vector of nodal values for each component

    no = ceil(sqrt(no));
    [X,Y,Wx,Wy] = TRIQUAD(no,V);    
    X = reshape(X,no^2,1);
    Y = reshape(Y,no^2,1);

    nouts  = nforces+nderivs;
    FJ     = cell(nouts,1);
    Felem  = cell(nforces,1);
    Jelem  = cell(nderivs,1);
    [FJ{:}] = func([X Y]);
    Nv     = TRI2D_SF([X Y]);
    
    for e=1:nforces
        Felem{e} = zeros(3,1);
        for n=1:3
            Felem{e}(n) = Wx'*(reshape(FJ{e}.*Nv(:,n),no,no))*Wy;
        end
    end
    ek = 1;
    for e=nforces+(1:nderivs)
        Jelem{ek} = zeros(3,3);
        for n=1:3
            for m=1:3
                Jelem{ek}(n,m) = Wx'*(reshape(Nv(:,n).*FJ{e}.*Nv(:,m),no,no))*Wy;
            end
        end
        ek = ek+1;
    end
end
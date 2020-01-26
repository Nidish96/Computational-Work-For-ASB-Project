function [Felem, Jelem] = QUAD2D_ELINT(func,V,no,nforces,nderivs)
%QUAD2D_ELINT Integrate given function over whole domain and report
%nodal components
% USAGE:
%	[Fnodal] = QUAD2D_ELINT(func,V,no,nforces,nderivs);
% INPUTS:
%   func	: (1x1)-valued vectorized function (in global CS)
%       	fval = func([x y]);
%       	  function could have multiple output scalars - see
%       	  nforces,nderivs
%   V 		: 4x2 vertex coordinates (global CS)
%   no		: (int) number of points for GQL - forced to be a
%   		  perfect square    
%   nforces	: Number of force output arguments of function
%   		  (first)
%   nderivs	: Number of force derivative output arguments of
%   		  function (last) 
% OUTPUTS:
%   Felem	: nforcesx1 cell of 4x1 vector of nodal values for
%   		  each component 
%   Jelem	: nderivsx1 cell of 4x1 vector of nodal values for
%   		  each component     
    
    no = ceil(sqrt(no));
    [X,Y,W,Jd] = QUADQUAD(no,V);
    X = reshape(X, no^2, 1);
    Y = reshape(Y, no^2, 1);
    
    nouts 	= nforces+nderivs;
    FJ		= cell(nouts,1);
    Felem	= cell(nforces,1);
    Jelem	= cell(nderivs,1);
    Nv		= QUAD2D_SF([X Y]);
    [FJ{:}]	= func(Nv*V);
    
    for e=1:nforces
        Felem{e} = zeros(4,1);
        for n=1:4
            Felem{e}(n) = W'*(reshape(FJ{e}.*Nv(:,n),no,no).*Jd)*W;
        end
    end
    ek = 1;
    for e=nforces+(1:nderivs)        
        Jelem{ek} = zeros(4,4);
        for n=1:4
            for m=1:4
                Jelem{ek}(n,m) = W'*(reshape(Nv(:,n).*FJ{e}.*Nv(:,m),no,no).*Jd)*W;
            end            
        end
        ek = ek+1;
    end
end
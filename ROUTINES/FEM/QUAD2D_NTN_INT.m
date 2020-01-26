function [NTNi] = QUAD2D_NTN_INT(V)
%QUAD2D_NTN_INT Returns the integral of the outer product of the
%shape functions vector against each other for a 2D quad defined by
%given node coordinates
% USAGE:
%	[NTNi] = QUAD2D_NTN_INT(V);
% INPUTS:
%   V		: 4x2 ordered coordinates of vertices (global CS)
% OUTPUTS:
%   NTNi	: 4x4 integrated matrix

% 2 GQL Points and weights
    w = [1.0 1.0];
    x = [-1.0/sqrt(3); 1.0/sqrt(3);];
    
    [xx,yy] = meshgrid(x,x);
    Ps = [reshape(xx,4,1) reshape(yy,4,1)];
    
    Ns = QUAD2D_SF(Ps);
    Js = QUAD2D_JACMAT(V,Ps);
    
    Jd = diag([det(Js(1:2,:)), det(Js(3:4,:)), det(Js(5:6,:)), ...
               det(Js(7:8,:))]); % All weights are unity, hence not
                                 % applied explicitly here
    if ~isempty(find(Jd<0, 1))
        error('Bad Jacobian - Quitting!');
    end
    NTNi = Ns'*Jd*Ns;
end
function [NTNi] = TRI2D_NTN_INT(V)
%TRI2D_NTN_INT Returns the integral of the outer product of the
%shape functions vector against each other for a 2D tri defined by
%given node coordinates
% USAGE:
%	[NTNi] = TRI2D_NTN_INT(V);
% INPUTS:
%   V 		: 3x2 ordered coordinates of vertices (global CS)
% OUTPUTS:
%   NTNi	: 3x3 integrated matrix
    
    % 3-point GQL
    w = [1.0; 1.0; 1.0]/3.0;
    Ps = [1.0/6 1.0/6;
          2.0/3 1.0/6;
          1.0/6 2.0/3];
    
    Ns = TRI2D_SF(Ps);
    Js = TRI2D_JACMAT(V,Ps);
    
    Jd = diag([det(Js(1:2,:))*w(1), det(Js(3:4,:))*w(2), det(Js(5:6,:))*w(3)]);
    NTNi = Ns'*Jd*Ns/2; % Check why this factor of 2 is required
end
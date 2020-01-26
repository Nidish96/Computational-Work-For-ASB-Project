function [Q12, Q21, N1TN1, N1TN2, N2TN2, N2TN1] = MESHTFM(MESH1, MESH2, No)
%MESHTFM returns the transformation matrix used to transfrom a DoF vector
%from MESH1 to MESH2 by minimizing the squared residual.
%
%  USAGE:
%       [Q, N1TN1, N2TN2, N1TN2, N2TN1] = MESHTFM(MESH1, MESH2, No);
%  INPUTS:
%   MESH1   : Mesh structure 1
%   MESH2   : Mesh structure 1
%   No      : Number of quadrature points
%  OUTPUTS:
%   Q12     : Transformation Matrix N1TN1\N1TN2; u1 = Q12.u2
%   Q21     : Transformation Matrix N1TN1\N1TN2; u2 = Q21.u1
%   N1TN1   : Integral of N1'*N1 over MESH1
%   N2TN2   : Integral of N2'*N2 over MESH2
%   N1TN2   : Integral of N1'*N2 over MESH2 
%   N2TN1   : Integral of N2'*N1 over MESH1

    % Operations on MESH 1
    [Q1, T1] = ZTE_ND2QP(MESH1, No);
	QP1 = Q1*MESH1.Nds;  % Quadrature points in MESH 1
    QP1_2 = MESH2PTS(MESH2, QP1);  % Interpolate MESH 2 onto QP 1
    
    N1TN1 = T1*Q1;
    N1TN2 = T1*QP1_2;
    
    % Final Transformation Matrix - FROM MESH 2 TO MESH 1 (all integrals on MESH 1)
    Q12 = N1TN1\N1TN2;   % u1 = Q12.u2 
    
    if nargout>1
        % Operations on MESH2
        [Q2, T2] = ZTE_ND2QP(MESH2, No);
        QP2 = Q2*MESH2.Nds;  % Quadrature points in MESH 2
        QP2_1 = MESH2PTS(MESH1, QP2);   % Interpolate MESH 1 onto QP 2
    
        N2TN2 = T2*Q2;
        N2TN1 = T2*QP2_1;
        
        % Final Transformation Matrix - FROM MESH 1 TO MESH 2 (all integrals on MESH 2)
        Q21 = N2TN2\N2TN1;   % u2 = Q21.u1
    end
end
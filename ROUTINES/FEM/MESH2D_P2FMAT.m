function [MAT] = MESH2D_P2FMAT(MESH)
%MESH2D_P2FMAT Returns the matrix transforming bi-linear
%distribution of pressure (collocated at nodes) to bi-linear
%distribution of forces. It effectively transforms nodal pressures
%into nodal forces.
% USAGE:
%	[MAT] = MESH2D_P2FMAT(MESH);
% INPUTS:
%   MESH	: Structure containing Mesh data
%       Nn      : Number of interface nodes
%       Nds     : (Nnx2) [x y] coordinates of interface nodes
%       Ne      : Number of interface elements
%       Ne_Tri	: Number of interface triangular elements
%       Ne_Quad	: Number of interface quadrilateral elements
%       Tri     : (Ne_Trix4) [eid n1 n2 n3]
%       Quad	: (Ne_Quadx4) [eid n1 n2 n3 n4]    
% OUTPUTS:
%   MAT		: (NnxNn) Transformation Matrix
    
    MAT = zeros(MESH.Nn,MESH.Nn);
    
    % Triangular Elements
    for e=MESH.Ne_Quad+(1:MESH.Ne_Tri)
        ek	= find(MESH.Tri(:,1)==e);
        nds 	= MESH.Tri(ek,2:4);
        V	= MESH.Nds(nds,:);
        MAT(nds,nds)	= MAT(nds,nds) + TRI2D_NTN_INT(V);
    end
    % Quadrilateral Elements
    for e=1:MESH.Ne_Quad
        ek 	= find(MESH.Quad(:,1)==e);
        nds	= MESH.Quad(ek,2:5);
        V	= MESH.Nds(nds,:);
        MAT(nds,nds)	= MAT(nds,nds) + QUAD2D_NTN_INT(V);
    end
end
function [P, Nums, NTNmat, NTGmat, GTGmat] = CONSPATCHMAT(Nds, Tri, Quad, XYvn)
%CONSPATCHMAT Returns the transformation from virtual DOFs (6)
%to original DOFs (3 each), consistent over a given patch defined
%by elements.
% USAGE:
%	[P, Nums] = CONSPATCHMAT(Nds, Tri, Quad, XYvn, no);
% INPUTS:
%   Nds		: Nnx2 locations of nodes
%   Tri		: Ntx4 triangular elements <eid, n1, n2, n3>
%   Quad	: Nqx4 quadrilateral elements <eid, n1, n2, n3, n4>
%   XYvn 	: 1x2 location of virtual node
%   no		: 1x1 Number of quadrature points in each direction
%   		  for integration
% OUTPUTS:
%   P		: Nn*3x6 transformation matrix
%   Nums    : Nn*3x1 "occurrence" vector
    Nn = size(Nds, 1);
    NTNmat = sparse(3*Nn, 3*Nn);
    NTGmat = sparse(3*Nn, 6);
    GTGmat = zeros(6, 6);
    
    Nt = size(Tri, 1);
    % Triangular Elements
    for e=1:Nt
        V = Nds(Tri(e, 2:end), :);
        is = reshape((Tri(e, 2:end)-1)*3 + (1:3)', 9, 1);
        NTNmat(is, is) = NTNmat(is, is) + ...
            sparse(kron(TRI2D_NTN_INT(V), eye(3)));
        
        XYnodal = TRI2D_ELINT(@(xy) deal(ones(size(xy, 1), 1), xy(:, 1)-XYvn(1), xy(:, 2)-XYvn(2)), V, 2, 3, 0);
        NTGmat(is(1:3:end), [1 6]) 	 = NTGmat(is(1:3:end), [1 6]) + [XYnodal{1} -XYnodal{3}];
        NTGmat(is(2:3:end), [2 6]) 	 = NTGmat(is(2:3:end), [2 6]) + [XYnodal{1} XYnodal{2}];
        NTGmat(is(3:3:end), [3 4 5]) = NTGmat(is(3:3:end), [3 4 5]) + [XYnodal{1} XYnodal{3} -XYnodal{2}];
    end
    
    Nq = size(Quad, 1);
    tm1 = zeros(3);
    tm2 = zeros(3);
    % Quadrilateral Elements
    for e=1:Nq
        V  = Nds(Quad(e, 2:end), :);
        is = reshape((Quad(e, 2:end)-1)*3 + (1:3)', 12, 1);
        NTNmat(is, is) = NTNmat(is, is) + ...
            sparse(kron(QUAD2D_NTN_INT(V), eye(3)));
        
        XYnodal = QUAD2D_ELINT(@(xy) deal(ones(size(xy,1),1), xy(:, 1)-XYvn(1), xy(:, 2)-XYvn(2)), V, 2, 3, 0);
        NTGmat(is(1:3:end), [1 6]) 	 = NTGmat(is(1:3:end), [1 6]) + [XYnodal{1} -XYnodal{3}];
        NTGmat(is(2:3:end), [2 6]) 	 = NTGmat(is(2:3:end), [2 6]) + [XYnodal{1} XYnodal{2}];
        NTGmat(is(3:3:end), [3 4 5]) = NTGmat(is(3:3:end), [3 4 5]) + [XYnodal{1} XYnodal{3} -XYnodal{2}];
        
        XY2nodal = QUAD2D_ELINT(@(xy) deal((xy(:, 1)-XYvn(1)).^2, (xy(:, 2)-XYvn(2)).^2, (xy(:, 1)-XYvn(1)).*(xy(:, 2)-XYvn(2))), ...
            V, 2, 3, 0);
        tm1 = [0 0 -sum(XYnodal{3}); 0 0 sum(XYnodal{2}); sum(XYnodal{3}) -sum(XYnodal{2}) 0];
        tm2 = [sum(XY2nodal{2}), -sum(XY2nodal{3}), 0;... 
            -sum(XY2nodal{3}), sum(XY2nodal{1}), 0;...
            0, 0, sum(XY2nodal{1}+XY2nodal{2})];
        GTGmat = GTGmat + [sum(XYnodal{1})*eye(3) tm1; tm1' tm2];
    end
    
    % Patch-Solving
    b0 = find(full(max(abs(NTGmat), [], 2))<eps);
    b1 = setdiff(1:size(NTGmat, 1), b0);
    P = zeros(Nn*3, 6);
    P(b1, :) = NTNmat(b1, b1)\NTGmat(b1, :);
    
    Nums = zeros(Nn*3, 1);
    Nums(b1) = 1;
end
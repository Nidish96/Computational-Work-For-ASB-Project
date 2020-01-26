function [A] = FEMGRAPHADJ(Nds,Quad_Els,Tri_Els)
%FEMGRAPHADJ Returns the adjacency matrix of the 2D mesh by having
%the element edges as node paths
% USAGE:
%	[A] = FEMGRAPHADJ(Nds,Quad_Els,Tri_Els);
% INPUTS:
%   Nds		: Nnx2 node locations
%   Tri_Els	: Ntx4 ordered list of nodes in Tri Elements
%		  [Elid n1 n2 n3]
%   Quad_Els	: Nqx5 ordered list of nodes in Quad Elements
%		  [Elid n1 n2 n3 n4]
% OUTPUTS:
%   A		: NnxNn Adjacency Matrix
    
    Nn = size(Nds,1);
    Nt = size(Tri_Els,1);
    Nq = size(Quad_Els,1);

    A = zeros(Nn,Nn);
    % Filling out Quads
    V = zeros(4,2);
    Vd = zeros(4,1);
    for e=1:Nq
        V = Nds(Quad_Els(e,2:5),:);
        Vd = sqrt(sum((diff(V([1:end 1],:))).^2,2));
        
        for k=1:4
            A(Quad_Els(e,1+k),Quad_Els(e,1+mod(k,4)+1)) = Vd(k);
            A(Quad_Els(e,1+mod(k,4)+1),Quad_Els(e,1+k)) = Vd(k);            
        end
    end
    % Filling out Tets
    V = zeros(3,2);
    Vd = zeros(3,1);
    for e=1:Nt
        V = Nds(Tri_Els(e,2:4),:);
        Vd = sqrt(sum((diff(V([1:end 1],:))).^2,2));
        
        for k=1:3
            A(Tri_Els(e,1+k),Tri_Els(e,1+mod(k,3)+1)) = Vd(k);
            A(Tri_Els(e,1+mod(k,3)+1),Tri_Els(e,1+k)) = Vd(k);
        end
    end
end
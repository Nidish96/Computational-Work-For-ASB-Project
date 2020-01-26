function [nid, eid] = MESH1_TO_MESH2(MESH1, MESH2)
%MESH1_TO_MESH2 Returns the node & elemenet ids of MESH1 corresponding to
%the nodes & elements ordered as per MESH2, such that
%MESH1.Nds(nid,:)=MESH2.Nds, and eid is the list of element ids in MESH1
%such that the elements are the same as MESH2.
%
% USAGE:
%   [nid, eid] = MESH1_TO_MESH2(MESH1, MESH2);
% INPUTS:
%   MESH1, MESH2   : Mesh 1 & 2 MESH Structure with elements
%       Nn      : Number of interface nodes
%       Nds     : (Nnx2) [x y] coordinates of interface nodes
%       Ne      : Number of interface elements
%       Ne_Tri	: Number of interface triangular elements
%       Ne_Quad	: Number of interface quadrilateral elements
%       Tri     : (Ne_Trix4) [eid n1 n2 n3]
%       Quad	: (Ne_Quadx4) [eid n1 n2 n3 n4]
%  OUTPUTS:
%   nid         : Node ids as per MESH1 so that they're along MESH2's nodes
%   eid         : Element ids as pert MESH2 so that they're along MESH2's
%                   elements
    % Normalizing both the meshes from 0 to 1
    Nd1 = MESH1.Nds;
    Nd1 = (Nd1 - min(Nd1))./range(Nd1);
    
    Nd2 = MESH2.Nds;
    Nd2 = (Nd2 - min(Nd2))./range(Nd2);
    
    % Nodes
    nid = zeros(size(Nd2, 1), 1);
    tmp = Nd1;
    for i=1:length(nid)
        [~, si] = sort(vecnorm(Nd1-Nd2(i,:),2,2));
        nid(i) = si(1);
    end
    
    % Elements sorted using centroids
    MESH1.Nds = Nd1; 
    MESH2.Nds = Nd2;
    
    [Q1, ~] = ZTE_ND2QP(MESH1, 1);
    [Q2, ~] = ZTE_ND2QP(MESH2, 1);
    
    Ec1 = Q1*Nd1;
    Ec2 = Q2*Nd2;
    
    eid = zeros(size(Ec2, 1),1);
    for i=1:length(eid)
        [~, si] = sort(vecnorm(Ec1-Ec2(i,:),2,2));
        eid(i) = si(1);
    end
end
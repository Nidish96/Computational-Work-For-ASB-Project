function [Q,Qx,Qy] = MESH2PTS(MESH, Pts)
%MESH2PTS returns a matrix that represents an interpolation from the nodes
%of the provided mesh to the required points. Brute-force searches employed
%here.
%
% USAGE:
%   Q = MESH2PTS(MESH, Pts);  % Verify: Pts = Q*MESH.Nds;
% INPUTS:
%   MESH        : Struct containing mesh information
%       Nn      : Number of nodes
%       Nds     : (Nnx2) [x y] coordinates of nodes
%       Ne      : Number of elements
%       Ne_Tri	: Number of triangular elements
%       Ne_Quad	: Number of quadrilateral elements
%       Tri     : (Ne_Trix4) [eid n1 n2 n3]
%       Quad	: (Ne_Quadx4) [eid n1 n2 n3 n4]
%   Pts         : (Npx2) matrix of points (x,y)
% OUTPUS:
%   Q           : NpxNn
    Pts = reshape(Pts, [], 2);
    Np = length(Pts);
    
    Tb = zeros(Np, MESH.Ne_Tri);  % Triangle Adjecency
    parfor e=1:MESH.Ne_Tri
%         Tb(:,e) = inpolygon(Pts(:,1), Pts(:,2), ...
%             MESH.Nds(MESH.Tri(e,[2:end 2]),1), MESH.Nds(MESH.Tri(e,[2:end 2]),2));
        
        Tb(:,e) = inpolygon(Pts(:,1), Pts(:,2), ...
            MESH.Nds(MESH.Tri(e,2:end),1), MESH.Nds(MESH.Tri(e,2:end),2));        
    end
    
    Qb = zeros(Np, MESH.Ne_Quad);  % Quad Adjecency
    parfor e=1:MESH.Ne_Quad
%         Qb(:,e) = inpolygon(Pts(:,1), Pts(:,2), ...
%             MESH.Nds(MESH.Quad(e,[2:end 2]),1), MESH.Nds(MESH.Quad(e,[2:end 2]),2));

        Qb(:,e) = inpolygon(Pts(:,1), Pts(:,2), ...
            MESH.Nds(MESH.Quad(e,2:end),1), MESH.Nds(MESH.Quad(e,2:end),2));
    end
    
    % Fill up matrix
    Q = zeros(Np, MESH.Nn);
    Qx = zeros(Np, MESH.Nn);
    Qy = zeros(Np, MESH.Nn);
    for p=1:Np
        % If in Triangle
        fT = find(Tb(p,:));
        if ~isempty(fT)
%             fT = find(MESH.Tri(:,1)==fT);
            for fti=1:length(fT)
                Q(p, MESH.Tri(fT(fti),2:end)) = Q(p, MESH.Tri(fT(fti),2:end))+TRI2D_SF_LOC(Pts(p,:), MESH.Nds(MESH.Tri(fT(fti),2:end),:));
            
                % Derivatives
                dxdy = TRI2D_SD_LOC(Pts(p,:), MESH.Nds(MESH.Tri(fT(fti),2:end),:));
                Qx(p, MESH.Tri(fT(fti),2:end)) = Qx(p, MESH.Tri(fT(fti),2:end)) + dxdy(1,:);
                Qy(p, MESH.Tri(fT(fti),2:end)) = Qx(p, MESH.Tri(fT(fti),2:end)) + dxdy(2,:);
            end
        end            
        
        % If in Quad
        fQ = find(Qb(p,:));
        if ~isempty(fQ)
%             fQ = find(MESH.Quad(:,1)==fQ(1));
%             QUAD2D_MAP2NAT(Pts(p,:), MESH.Nds(MESH.Quad(fQ, 2:end),:)); % Mapped

            for fqi=1:length(fQ)
                Q(p, MESH.Quad(fQ(fqi),2:end)) = Q(p, MESH.Quad(fQ(fqi),2:end)) + QUAD2D_SF(QUAD2D_MAP2NAT(Pts(p,:), MESH.Nds(MESH.Quad(fQ(fqi), 2:end),:)));
                
                % Derivatives
                dxideta = QUAD2D_SD(QUAD2D_MAP2NAT(Pts(p,:), MESH.Nds(MESH.Quad(fQ(fqi), 2:end),:)));
                [~,Ji] = QUAD2D_JACMAT(MESH.Nds(MESH.Quad(fQ(fqi), 2:end),:), Pts(p,:));
                Qx(p, MESH.Quad(fQ(fqi),2:end)) = Qx(p, MESH.Quad(fQ(fqi),2:end)) + Ji(:,1)'*dxideta;
                Qy(p, MESH.Quad(fQ(fqi),2:end)) = Qy(p, MESH.Quad(fQ(fqi),2:end)) + Ji(:,2)'*dxideta;
            end
        end
        
        if length(fQ)+length(fT)>0
            Q(p, :) = Q(p, :)/(length(fQ)+length(fT));
            Qx(p, :) = Qx(p, :)/(length(fQ)+length(fT));
            Qy(p, :) = Qy(p, :)/(length(fQ)+length(fT));
        end
    end
    
    % Find closest node for all points that don't seem to be within an
    % element
    nf = find(~sum(Q~=0,2));
    dist = sqrt((Pts(:,1)-MESH.Nds(:,1)').^2 + (Pts(:,2)-MESH.Nds(:,2)').^2);
    for p=1:length(nf)
        [~,mp] = min(dist(nf(p),:));
        Q(nf(p), mp) = 1;
    end
end
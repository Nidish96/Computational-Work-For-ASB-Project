function [] = SHOW3D(Nds,Tri_Els,Quad_Els,Ux,Uy,Uz,sc,P)
%SHOW3D Plots the elements in 3D with given x,y & z displacements
% USAGE:
%	[] = SHOW3D(Nds,Tri_Els,Quad_Els,Ux,Uy,Uz,sc);
% INPUTS:
%   Nds,Tri_Els,Quad_Els,Ux,Uy,Uz,sc
% OUTPUTS:
%   

    Nn  = size(Nds, 1);
    Nt 	= size(Tri_Els,1);
    Nq	= size(Quad_Els,1);
    Ne = Nt + Nq;
    
    % Triangles
    for e=1:Nt
        nds	= Tri_Els(e,2:4);
        V	= Nds(nds,:);
%         fill3(V(:,1)+Ux(nds)*sc,V(:,2)+Uy(nds)*sc,Uz(nds)*sc,P);
        
        if size(P,1)==Nn
            fill3(V(:,1)+Ux(nds)*sc,V(:,2)+Uy(nds)*sc,Uz(nds)*sc,P(Tri_Els(e,2:end),:));
        elseif size(P,1)==Ne
            fill3(V(:,1)+Ux(nds)*sc,V(:,2)+Uy(nds)*sc,Uz(nds)*sc,P(Tri_Els(e,1),:))
        else
            fill3(V(:,1)+Ux(nds)*sc,V(:,2)+Uy(nds)*sc,Uz(nds)*sc,P(:,1));
        end
        hold on;
    end
    % Quadrilaterals
    for e=1:Nq
        nds	= Quad_Els(e,2:5);
        V	= Nds(nds,:);
%         fill3(V(:,1)+Ux(nds)*sc,V(:,2)+Uy(nds)*sc,Uz(nds)*sc,P);
        
        if size(P,1)==Nn
            fill3(V(:,1)+Ux(nds)*sc,V(:,2)+Uy(nds)*sc,Uz(nds)*sc,P(Quad_Els(e,2:end),:));
        elseif size(P,1)==Ne
            fill3(V(:,1)+Ux(nds)*sc,V(:,2)+Uy(nds)*sc,Uz(nds)*sc,P(Quad_Els(e,1),:))
        else
            fill3(V(:,1)+Ux(nds)*sc,V(:,2)+Uy(nds)*sc,Uz(nds)*sc,P(1,:));
        end
        hold on;
    end
end
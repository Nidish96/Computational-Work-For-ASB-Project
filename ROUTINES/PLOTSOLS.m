function [] = PLOTSOLS(Nds,Tri_Els,Quad_Els,Sols, Titles, varargin)
%PLOTSOLS Plots provided functions in provided cell format
% USAGE:
%	[] = PLOTSOLS(Nds,Tri_Els,Quad_Els,Sols);
% INPUTS:
%   Nds		: (Nnx2)
%   Tri_Els	: (Ntx2) 
%   Quad_Els	: (Nqx2) 
%   Sols 	: {rows x cols}
% OUTPUTS:
%   
    if nargin < 6
        cpos = 'southoutside';
        cscl = 'linear';
        outl = 0;
        fsz = 9;
    elseif nargin == 6
        cpos = varargin{1};
        cscl = 'linear';
        outl = 0;
        fsz = 9;
    elseif nargin == 7
        cpos = varargin{1};
        cscl = varargin{2};
        outl = 0;
        fsz = 9;
    elseif nargin == 8
        cpos = varargin{1};
        cscl = varargin{2};
        outl = varargin{3};
        fsz = 9;
    elseif nargin == 9
        cpos = varargin{1};
        cscl = varargin{2};
        outl = varargin{3};
        fsz = varargin{4};
    end

    rows = size(Sols,1);
    cols = size(Sols,2);
    Ntot = rows*cols;
    
    Ntri = size(Tri_Els,1);
    Nquad = size(Quad_Els,1);
    for i=1:rows
        for j=1:cols
            subplot(rows,cols, (i-1)*cols+j)
            
            % Triangles
            for e=1:Ntri
                nds = Tri_Els(e,2:4);
                V   = Nds(nds,:);
                if length(Sols{i,j})==size(Nds,1)
                    fill(V(:,1),V(:,2),Sols{i,j}(nds)); hold on
                elseif length(Sols{i,j})==(Ntri+Nquad)
                    fill(V(:,1),V(:,2),Sols{i,j}(e)); hold on
                end
                if outl==1
                    plot(V([1:end 1],1),V([1:end 1],2), 'k-')
                end
            end
            % Quadrilaterals
            for e=1:Nquad
                nds = Quad_Els(e,2:5);
                V   = Nds(nds,:);
                if length(Sols{i,j})==size(Nds,1)
                    fill(V(:,1),V(:,2),Sols{i,j}(nds)); hold on
                elseif length(Sols{i,j})==(Ntri+Nquad)
                    fill(V(:,1),V(:,2),Sols{i,j}(e)); hold on
                end
                if outl==1
                    plot(V([1:end 1],1),V([1:end 1],2), 'k-')
                end
            end
            
%             axis equal
            axis off
            cb = colorbar(cpos);
            title(Titles{i,j}, 'FontWeight', 'Normal', 'interpreter', 'latex')
            set(gca, 'ColorScale', cscl, 'FontSize', fsz);
            if strcmp(cscl, 'log')
%                 caxis(10.^[floor(log10(min(Sols{i,j}))) ceil(log10(max(Sols{i,j})))])
%                 set(cb, 'Ticks', logspace(floor(log10(min(Sols{i,j}))), ceil(log10(max(Sols{i,j}))), 4)) 
                set(cb, 'Ticks', logspace(-20, 20, 41))
            end
%             set(gca, 'FontSize', fsz);
        end
    end    
end
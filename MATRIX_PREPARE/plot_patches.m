clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/FEM/')
addpath('../ROUTINES/export_fig/')

conf = 'CD1'; % 'C1' 'C2' 'C3' 'C4' 'C4a' 'CD1'

Confs = {'C1', 'C2', 'C3', 'C4', 'C4a', 'CD1'};
for ci=1:length(Confs)
    conf = Confs{ci};
%% Load
load(sprintf('./MAT_NULLRED_%s.mat', conf), 'Pels', 'MESH', 'Npatches')

%% Plot
colos = (1:Npatches)';

figure(1)
clf()
for i=1:Npatches
    SHOW2DMESH(MESH.Nds, MESH.Tri, MESH.Quad(Pels{i},:), colos(i,:), -1, -100, 0, MESH.Ne)
end
axis equal
axis off

export_fig(sprintf('./FIGS/PATCHES_%s.png', conf), '-dpng','-r1200')
end
clc
clear all 
addpath('../ROUTINES')
addpath('../ROUTINES/FEM')

conf = '4a'; 
mfname = sprintf('../MATRIX_PREPARE/MAT_NULLRED_%s.mat',conf);
load(mfname, 'M', 'K', 'L', 'R', 'Fv', 'Th', 'LamT', 'cnum', 'MESH', 'PatchAreas', 'Npatches', 'dofred', 'GTG')

%% Access Matrices and Functions
Qxyn = L(reshape(((1:Npatches)-1)*6+(1:3)', Npatches*3, 1), :);
Txyn = LamT(:, reshape(((1:Npatches)-1)*6+(1:3)', Npatches*3, 1));

%% Parameters
kt = 2.5e8;
kn = 4e8;
Kmat = diag(kron(ones(1, Npatches), [kt kt kn]));

ct = 1e5;
cn = 1e5;
Cmat = diag(kron(ones(1, Npatches), [ct ct cn]));

K0 = K + Txyn*Kmat*Qxyn;
C0 = Txyn*Cmat*Qxyn;

[V, Ws] = eigs(K0, M, 10, 'SM');
[Ws, si] = sort(sqrt(diag(Ws)));
V = V(:, si);
V = V./sqrt(diag(V'*M*V)');

%% Linear Modal Properties
Zs = diag(V'*C0*V)./(2*Ws);
Ws = Ws/2/pi;

disp(Ws')
disp(Zs')

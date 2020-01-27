clc
clear all 
addpath('../ROUTINES')
addpath('../ROUTINES/FEM')

conf = '1'; 
mfname = sprintf('../MATRIX_PREPARE/MAT_NULLRED_C%s.mat',conf);
load(mfname, 'M', 'K', 'L', 'R', 'Fv', 'Th', 'LamT', 'cnum', 'MESH', 'PatchAreas', 'Npatches', 'dofred', 'GTG')

%% Access Matrices and Functions
Qxyn = L(reshape(((1:Npatches)-1)*6+(1:3)', Npatches*3, 1), :);
Txyn = LamT(:, reshape(((1:Npatches)-1)*6+(1:3)', Npatches*3, 1));

Txyn = Qxyn';

opts.dK = cell(4, 1);
opts.dK{1} = Txyn*diag(kron(ones(1, Npatches), [1 1 0]))*Qxyn;
opts.dK{2} = Txyn*diag(kron(ones(1, Npatches), [0 0 1]))*Qxyn;

opts.dC = cell(4, 1);
opts.dC{3} = Txyn*diag(kron(ones(1, Npatches), [1 1 0]))*Qxyn;
opts.dC{4} = Txyn*diag(kron(ones(1, Npatches), [0 0 1]))*Qxyn;

opts.lspci = [];

%% Parameters
kt = 2e8;
kn = 4e8;
ct = 1e5;
cn = 1e5;

pars = [kt kn ct cn];
pars(opts.lspci) = log10(pars(opts.lspci));

%% Function Call and Derivatives
h = 1e-4;

for hi = 1:4
    hv = zeros(size(pars));
    hv(hi) = 1;


    [s, dsdpar] = QEP_ROOTS(pars, K, M*0, M, 1, opts);
    sp = QEP_ROOTS(pars+hv*h, K, M*0, M, 1, opts);
    sm = QEP_ROOTS(pars-hv*h, K, M*0, M, 1, opts);
    fprintf('Derivative %d:', hi)
    disp([(sp-sm)/(2*h) dsdpar(hi)])
end

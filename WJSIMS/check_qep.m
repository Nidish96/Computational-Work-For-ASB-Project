clc
clear all 
addpath('../ROUTINES')
addpath('../ROUTINES/FEM')

conf = '1'; 
mfname = sprintf('../MATRIX_PREPARE/MAT_NULLRED_C%s.mat',conf);
load(mfname, 'M', 'K', 'L', 'R', 'Fv', 'Th', 'LamT', 'cnum', 'MESH', 'PatchAreas', 'Npatches', 'dofred', 'GTG')

%% Experimental Data
exx = load(sprintf('../EXP_DATA/ASB_C%s_20Nm_BB.mat', conf), 'W', 'Z');

%% Access Matrices and Functions
Qxyn = L(reshape(((1:Npatches)-1)*6+(1:3)', Npatches*3, 1), :);
Txyn = LamT(:, reshape(((1:Npatches)-1)*6+(1:3)', Npatches*3, 1));

opts.dK = cell(4, 1);
opts.dK{1} = Txyn*diag(kron(ones(1, Npatches), [1 1 0]))*Qxyn;
opts.dK{2} = Txyn*diag(kron(ones(1, Npatches), [0 0 1]))*Qxyn;

opts.dC = cell(4, 1);
opts.dC{3} = Txyn*diag(kron(ones(1, Npatches), [1 1 0]))*Qxyn;
opts.dC{4} = Txyn*diag(kron(ones(1, Npatches), [0 0 1]))*Qxyn;

opts.lspci = [1:4];

opts.mode = 'standard'; % {'standard', 'polar', 'sods'}
opts.qtty = 'sq_dev'; % {'abs_dev', 'sq_dev', 'value'}

opts.lspco = [1:2];

opts.expdat = [2*pi*mean(exx.W); mean(exx.Z)];

%% Parameters
kt = 2e8;
kn = 4e8;
ct = 1e5;
cn = 1e5;

pars = [kt kn ct cn];
pars(opts.lspci) = log10(pars(opts.lspci));

%% Function Call and Derivative Check
h = 1e-7;

opts.mode = 'sods';
for hi = 1:4
    hv = zeros(size(pars));
    hv(hi) = 1;

    [s, dsdpar] = QEP_ROOTS(pars, K, M*0, M, 1, opts);
    sp = QEP_ROOTS(pars+hv*h, K, M*0, M, 1, opts);
    sm = QEP_ROOTS(pars-hv*h, K, M*0, M, 1, opts);
    fprintf('Derivative %d:\n', hi)
				% disp([(sp-sm)/(2*h) dsdpar(hi)])
    nders = [(sp-sm)/2 (sp-s) (s-sm)]/(hv(hi)*h);
    switch (opts.mode)
      case 'standard'
	fprintf('(%.4e,%.4e); (%.4e,%.4e); (%.4e,%.4e). (%.4e,%.4e)\n', ...
		abs(nders(2)), angle(nders(2)), ...
		abs(nders(1)), angle(nders(1)), ...
		abs(nders(3)), angle(nders(3)), ...
		abs(dsdpar(hi)), angle(dsdpar(hi)))
      case {'polar', 'sods'}
	fprintf('(%.4e,%.4e); (%.4e,%.4e); (%.4e,%.4e). (%.4e,%.4e)\n', ...
		nders(1, 1), nders(2, 1), ...
		nders(1, 2), nders(2, 2), ...
		nders(1, 3), nders(2, 3), ...
		dsdpar(1, hi), dsdpar(2, hi))
    end
end

%% Genetic Algorithm
gopts = optimoptions('gamultiobj', 'Display', 'iter', 'PlotFcn', ...
                     {@gaplotpareto, @gaplotrankhist}, 'MaxGenerations', ...
                     100, 'DistanceMeasureFcn', {@distancecrowding, ...
                    'genotype'}, 'PopulationSize', 100); %,
                                                         %'InitialPopulationMatrix',
                                                         %[PARS]);
lb = [-1, -1, -1, -1];
ub = [30, 30, 30, 30];
tic
[PARS_GA, ERRS_GA] = gamultiobj(@(prs) QEP_ROOTS(prs, K, M*0, M, 1, opts), length(pars), ...
                                [], [], [], [], lb, ub, [], gopts);
toc
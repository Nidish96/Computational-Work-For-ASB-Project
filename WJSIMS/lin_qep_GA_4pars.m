clc
clear all 
addpath('../ROUTINES')
addpath('../ROUTINES/FEM')

conf = '1';
Confs = {'1', '2', '3', '4', '4a', 'D1'};

for ic=1:length(Confs)
    conf = Confs{ic};
    
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

    opts.lspci = 1:4;

    opts.mode = 'sods'; % {'standard', 'polar', 'sods'}
    opts.qtty = 'sq_dev'; % {'abs_dev', 'sq_dev', 'value'}

    opts.lspco = 1:2;

    opts.expdat = [2*pi*mean(exx.W); mean(exx.Z)];

    %% Genetic Algorithm
    gopts = optimoptions('gamultiobj', 'Display', 'iter', 'PlotFcn', ...
                         {@gaplotpareto, @gaplotrankhist}, 'MaxGenerations', ...
                         200, 'DistanceMeasureFcn', {@distancecrowding, ...
                        'genotype'}, 'PopulationSize', 500, 'UseParallel', true); %,
                                                             %'InitialPopulationMatrix',
                                                             %[PARS]);
    lb = [-1, -1, -1, -1];
    ub = [30, 30, 30, 30];
    tic
    [PARS_GA, ERRS_GA] = gamultiobj(@(prs) QEP_ROOTS(prs, K, M*0, M, 1, opts), length(lb), ...
                                    [], [], [], [], lb, ub, [], gopts);
    toc

    save(sprintf('./DATS/ASB_C%s_LINQEP_GA_RES.mat', conf), 'PARS_GA', 'ERRS_GA');
end
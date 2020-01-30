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
%     opts.lspco = [];

    opts.expdat = [2*pi*mean(exx.W); mean(exx.Z)];

    %% Get Utopian point from GA runs
    load(sprintf('./DATS/ASB_C%s_LINQEP_GA_RES.mat', conf), 'PARS_GA', 'ERRS_GA');
    [Wmin, iw] = min(ERRS_GA(:,1));  PARMIN_W = PARS_GA(iw,:);
    [Zmin, iz] = min(ERRS_GA(:,2));  PARMIN_Z = PARS_GA(iz,:);
    
    %% Spherical Weighting
    sopts = struct('method', 'sphericalbi', 'gradients', true, ...
        'nobj', 2, 'npar', 2, 'nvar', length(opts.dK));
    sopts.rpt = [Wmin; Zmin]; % - abs([Wmin; Zmin])*0.1;
    sopts.rpt = [-25; -20];
%     sopts.rpt = [0; 0];
    
    THETAS_BI = linspace(0, pi/2, 50);
    ERRS_BI = zeros(length(THETAS_BI), 2);
    PARS_BI = zeros(length(THETAS_BI), sopts.nvar+1);
    EXITFLAG = zeros(length(THETAS_BI), 1);
    
%     pars0 = (PARMIN_W(:)+PARMIN_Z(:))/2;
    pars0 = PARMIN_W';
    lb = [-1, -1, -1, -1, 0];
    ub = [30, 30, 30, 30, +Inf];
%     pars0 = (lb(:)+ub(:))/2;
%     pars0 = mean(PARS_GA)';
    
    fopts = optimoptions('fmincon', 'SpecifyObjectiveGradient', true, ...
        'Display', 'iter', 'SpecifyConstraintGradient', false, ...
        'MaxIterations', 400, 'StepTolerance', 1e-20);
    fopts.Algorithm = 'interior-point';
    t0 = 6;
    
    for k=1:length(THETAS_BI)
        [PARS_BI(k,:), err, EXITFLAG(k)] = ...
            fmincon(@(prs) deal(prs(end), [zeros(1, length(prs)-1) 1]), ...
            [pars0; t0], [], [], [], [], lb, ub, ...
            @(prs) SCALARIZE(@(pars) QEP_ROOTS(pars, K, M*0, M, 1, opts), ...
            [prs; THETAS_BI(k)], sopts), fopts);
        ERRS_BI(k,:) = QEP_ROOTS(PARS_BI(k, 1:end-1), K, M*0, M, 1, opts);
        
        fprintf('Done %d/%d\n', k, length(THETAS_BI));
        
        t0 = PARS_BI(k, end);
    end
    
    clf(); plot(ERRS_BI(:,1), ERRS_BI(:,2), 'o', ERRS_GA(:,1), ERRS_GA(:,2), 'x', Wmin, Zmin, 'X')

    save(sprintf('./DATS/ASB_C%s_LINQEP_SBI_RES.mat', conf), 'PARS_GA', 'ERRS_GA');
end
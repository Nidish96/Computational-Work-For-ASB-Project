clc
clear all 
addpath('../ROUTINES')
addpath('../ROUTINES/WJSIMS')
addpath('../ROUTINES/WJSIMS/IWAN')
addpath('../ROUTINES/FEM')
addpath('../ROUTINES/IWAN')

conf = '1';
Confs = {'1', '2', '3', '4', '4a', 'D1'};

for ci=1:length(Confs)  % 3,4,6 not running
    conf = Confs{ci};
    mfname = sprintf('../MATRIX_PREPARE/MAT_NULLRED_C%s.mat',conf);
    load(mfname, 'M', 'K', 'L', 'R', 'Fv', 'Th', 'LamT', 'cnum', ...
         'MESH', 'PatchAreas', 'Npatches', 'dofred', 'GTG')
    Prestress = 11580;

    %% Access Matrices and Functions
    % LamT = LamT*GTG;
    Qxyn = L(reshape(((1:Npatches)-1)*6+(1:3)', Npatches*3, 1), :);
    Txyn = LamT(:, reshape(((1:Npatches)-1)*6+(1:3)', Npatches*3, 1));

    copt.lspci = [1 2 4 5 6];  % [Fs Kt Chi Bt Kp Kn]  (Chi not log-scale)

    copt.x.T{1} = [PatchAreas(:) zeros(Npatches, 5)];  % Fs
    copt.x.T{2} = [zeros(Npatches,1) PatchAreas(:) zeros(Npatches, 4)];  % Kt
    copt.x.T{3} = [zeros(Npatches,2) ones(Npatches,1) zeros(Npatches, 3)];  % Chi
    copt.x.T{4} = [zeros(Npatches,3) ones(Npatches,1) zeros(Npatches, 2)];  % Bt
    copt.x.T{5} = [zeros(Npatches,4) PatchAreas(:) zeros(Npatches, 1)];  % Kp

    copt.y.T{1} = [PatchAreas(:) zeros(Npatches, 5)];  % Fs
    copt.y.T{2} = [zeros(Npatches,1) PatchAreas(:) zeros(Npatches, 4)];  % Kt
    copt.y.T{3} = [zeros(Npatches,2) ones(Npatches,1) zeros(Npatches, 3)];  % Chi
    copt.y.T{4} = [zeros(Npatches,3) ones(Npatches,1) zeros(Npatches, 2)];  % Bt
    copt.y.T{5} = [zeros(Npatches,4) PatchAreas(:) zeros(Npatches, 1)];  % Kp

    copt.n.T{1} = [zeros(Npatches, 5) PatchAreas(:)];  % Fs

    CFUN = @(uxyn, pars) SPRING_IWAN4_MOD(uxyn, pars, copt);
    DFUN = @(uxyn, pars) SPRING_IWAN4_MOD_DISS(uxyn, pars, copt);

    %% Load Experimental Data
    exx = load(sprintf('../EXP_DATA/ASB_C%s_20Nm_BB.mat', conf), ...
               'Q', 'W', 'Z');
    Nx = 20; iNs = fix(linspace(1, length(exx.Q), Nx));
    expdat.Q = exx.Q(iNs);
    expdat.W = exx.W(iNs)*2*pi;
    expdat.Z = exx.Z(iNs);
    expdat.D = 2*pi*(expdat.Q.*expdat.W).^2.*expdat.Z;
    
    %% Get Utopian point from GA runs
    load(sprintf('./DATS/ASB_C%s_4IWAN_PSSTIFF_GA_RES.mat', conf), ...
         'PARS_GA', 'ERRS_GA')
    [Wmin, iw] = min(ERRS_GA(:,1));  PARMIN_W = PARS_GA(iw,:);
    [Zmin, iz] = min(ERRS_GA(:,2));  PARMIN_Z = PARS_GA(iz,:);
    
    % Fully Stuck Initial Guess with PARMIN_W
    lpars = PARMIN_W(:); lpars(copt.lspci) = 10.^(lpars(copt.lspci)); 
    Kst = Txyn*diag(reshape([copt.x.T{2}*lpars copt.y.T{2}*lpars ...
                        copt.n.T{1}*lpars]', Npatches*3, 1))*Qxyn;
    K0 = K + Kst;
    X0 = K0\(Prestress*Fv);

    %% Spherical Weighting
    opts = struct('method', 'sphericalbi', 'gradients', true, ...
                  'nobj', 2, 'npar', 2, 'nvar', size(copt.x.T{1},2));
    opts.rpt = [Wmin; Zmin];
    opt = optimoptions('fsolve', 'Display', 'off', ...
                       'SpecifyObjectiveGradient', true);
    % opts.rpt = [0; 0];
    Nqp = 50; 

    THETAS_BI = linspace(0, pi/2, 10);
    % THETAS_BI = [0 0.001 0.005 0.01 0.015 0.02 0.025 0.035 0.05 0.075 ...
    %              0.10 0.15 0.20 0.25 0.75 1.0]*pi/2;
    ERRS_BI = zeros(length(THETAS_BI), 2);
    PARS_BI = zeros(length(THETAS_BI), opts.nvar+1);
    EXITFLAG = zeros(length(THETAS_BI), 1);

    PARS_BI(1,:) = [PARMIN_W(:); 0.0];
    ERRS_BI(1,:) = WJMODEL_BBFUN(PARS_BI(1,1:end-1), 1, expdat, K, M, X0, ...
                                 Fv*Prestress, L, Txyn, Qxyn, CFUN, ...
                                 Npatches, Nqp, opt);

    pars0 = (PARMIN_W(:) + PARMIN_Z(:))/2;
    
    lb = [0, 0, -2, -6, 0, 3];
    ub = [30, 20, 0, 5, 22, 20];
    
    fopts = optimoptions('fmincon', 'SpecifyObjectiveGradient', true, ...
                         'Display', 'iter', 'SpecifyConstraintGradient', ...
                         true, 'MaxIterations', 400);
    t0 = 2;
    
    for k=1:length(THETAS_BI)
        [PARS_BI(k,:), err, EXITFLAG(k)] = fmincon(@(prs) deal(prs(end), ...
                                                          [zeros(1, ...
                                                          length(prs)-1) ...
                            1]), [pars0; t0], [], [], [], [], lb, ...
                                                   ub, @(prs) ...
                                                   SCALARIZE(@(pars) ...
                                                          WJMODEL_BBFUN_MASING(pars, ...
                                                          1, expdat, ...
                                                          K, M, X0, ...
                                                          Fv*Prestress, ...
                                                          L, Txyn, ...
                                                          Qxyn, CFUN, ...
                                                          DFUN, ...
                                                          Npatches, ...
                                                          opt), [prs; ...
                            THETAS_BI(k)], opts), fopts);
        ERRS_BI(k,:) = WJMODEL_BBFUN(PARS_BI(k,1:end-1), 1, expdat, K, M, ...
                                     X0, Fv*Prestress, L, Txyn, Qxyn, ...
                                     CFUN, Npatches, Nqp, opt);
        
        fprintf('Done %d/%d\n', k, length(THETAS_BI));
    end

    save(sprintf('./DATS/ASB_C%s_4IWAN_PSSTIFF_BI_RES.mat', conf), 'PARS_BI', 'ERRS_BI', ...
         'PARMIN_W', 'PARMIN_Z');

    % figure(10)
    % % clf()
    plot(ERRS_BI(:,1), ERRS_BI(:,2), '.-'); hold on
end

%% Backbone function
opt.Display = 'off';

pars = [5.5; 13.8; -0.95; -10; 11.6773; 12.35]
pars = PARMIN_Z;
Nqp = 50; 
tic
% [eobj, dedp, BB] = WJMODEL_BBFUN(pars, 1, expdat, K, M, X0, ...
%                                  Fv*Prestress, L, Txyn, Qxyn, ...
%                                  CFUN, Npatches, Nqp, opt);
[eobj_2, dedp_2, BB_2] = WJMODEL_BBFUN_MASING(pars, 1, expdat, ...
                                              K, M, X0, Fv*Prestress, ...
                                              L, Txyn, Qxyn, CFUN, DFUN, Npatches, opt);
toc

% figure(10) 
% % clf()
% subplot(2,1,1)
% semilogx(BB_2.Q, BB_2.W/(2*pi), '.-', expdat.Q, expdat.W/(2*pi), ...
%          'k-'); hold on
% subplot(2,1,2)
% loglog(BB_2.Q, BB_2.Z, '.-', expdat.Q, expdat.Z, 'k-'); hold on
% disp('Done')
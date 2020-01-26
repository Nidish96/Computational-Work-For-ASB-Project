clc
clear all 
addpath('../ROUTINES')
addpath('../ROUTINES/WJSIMS')
addpath('../ROUTINES/WJSIMS/IWAN')
addpath('../ROUTINES/FEM')
addpath('../ROUTINES/IWAN')

conf = '1';
Confs = {'1', '2', '3', '4', '4a', 'D1'};

for ci=1:length(Confs)
    conf = Confs{ci};
    mfname = sprintf('../MATRIX_PREPARE/MAT_NULLRED_%s.mat',conf);
    load(mfname, 'M', 'K', 'L', 'R', 'Fv', 'Th', 'LamT', 'cnum', 'MESH', 'PatchAreas', 'Npatches', 'dofred', 'GTG')
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

    %% Contact Model
    pars = [9; 20; -0.8; -4; 11.6773; 20];  % [Fs Kt Chi Bt Kn]
    lpars = pars; lpars(copt.lspci) = 10.^(lpars(copt.lspci)); 

    %% Fully Stuck Initial Guess For Prestress Analysis
    Kst = Txyn*diag(reshape([copt.x.T{2}*lpars copt.y.T{2}*lpars ...
                        copt.n.T{1}*lpars]', Npatches*3, 1))*Qxyn;
    K0 = K + Kst;
    X0 = K0\(Prestress*Fv);

    %% Prestress Analysis
    opt = optimoptions('fsolve', 'Display', 'iter', 'SpecifyObjectiveGradient', true);

    [Xstat, ~, eflg] = fsolve(@(X) WJRES_STR([X; 0.0], pars, K, Fv*0, ...
                                             Fv*Prestress, L, Txyn, ...
                                             Qxyn, CFUN, Npatches), ...
                              X0, opt);
    [~, dRstat, ~, dRdpstat] = WJRES_STR([Xstat; 0.0], pars, K, Fv*0, ...
                                         Fv*Prestress, L, Txyn, Qxyn, ...
                                         CFUN, Npatches);
    dXdpstat = -dRstat\dRdpstat;

    %% Modal Analysis
    [Vst, Wst] = eigs(dRstat, M, 10, 'SM');
    Wst = sqrt(diag(Wst));
    Vst = Vst./sqrt(diag(Vst'*M*Vst)');

    %% Load Experimental Data
    exx = load(sprintf('../EXP_DATA/ASB_C%s_20Nm_BB.mat', conf), ...
               'Q', 'W', 'Z');
    Nx = 20; iNs = fix(linspace(1, length(exx.Q), Nx));
    expdat.Q = exx.Q(iNs);
    expdat.W = exx.W(iNs)*2*pi;
    expdat.Z = exx.Z(iNs);
    expdat.D = 2*pi*(expdat.Q.*expdat.W).^2.*expdat.Z;

    %% Scalarized Optimization for utopian point
    opts = struct('method', 'weighted', 'gradients', true, 'nobj', ...
        2, 'npar', 1, 'nvar', length(pars));
    fopts = optimoptions('fmincon', 'SpecifyObjectiveGradient', true, ...
        'Display', 'iter');
    opt.Display = 'off';
    % gopts = optimoptions('ga', 'Display', 'iter', 'MaxGenerations', 20);
    lb = [3, 3, -1, -4, 0, 3];
    ub = [30, 20, 0, 4, 22, 20];

    % pars0 = [1, 2, -0.5, 3, 4]';
    pars0 = [5; 13.8; -0.5; -4; 11.6773; 12.35];

    wgt = 1.0;  % Frequency Alone
                % Masing/Diss function forms
    [PARMIN_W, Wmin, EXITFLAG_W] = fminunc(@(prs) SCALARIZE(@(pars) WJMODEL_BBFUN_MASING(pars, 1, expdat, K, M, X0, Fv*Prestress, L, Txyn, Qxyn, CFUN, DFUN, Npatches, opt), [prs; wgt], opts), ...
                                           pars0, fopts);  % unconstrained opt
    wgt = 0.0;  % Damping Alone
                %Masing/Diss function forms
    [PARMIN_Z, Zmin, EXITFLAG_Z] = fminunc(@(prs) SCALARIZE(@(pars) ...
                                                      WJMODEL_BBFUN_MASING(pars, 1, expdat, K, M, X0, Fv*Prestress, L, Txyn, Qxyn, CFUN, DFUN, Npatches, opt), [prs; wgt], opts), pars0, fopts); % unconstrained opt

    % Utopian Point
    opts.rpt = [Wmin; Zmin];

    %% Spherical Weighting
    opts = struct('method', 'sphericalwgts', 'gradients', true, ...
                  'nobj', 2, 'npar', 1, 'nvar', length(pars));  
    % opts.rpt = [Wmin; Zmin];
    opts.rpt = [0; 0];
    Nqp = 50; 

    THETAS_SW = linspace(0, pi/2, 10);
    % THETAS_SW = [0 0.001 0.005 0.01 0.015 0.02 0.025 0.035 0.05 0.075 ...
    %              0.10 0.15 0.20 0.25 0.75 1.0]*pi/2;
    ERRS_SW = zeros(length(THETAS_SW), 2);
    PARS_SW = zeros(length(THETAS_SW), length(pars));

    PARS_SW(1,:) = PARMIN_W;
    ERRS_SW(1,:) = WJMODEL_BBFUN(PARS_SW(1,:), 1, expdat, K, M, X0, ...
                                 Fv*Prestress, L, Txyn, Qxyn, CFUN, Npatches, Nqp, opt);

    pars0 = (PARMIN_W + PARMIN_Z)/2;
    
    fopts.Display = 'off';
    opt.Display = 'off';
    for k=1:length(THETAS_SW)
        [PARS_SW(k, :), err] = fminunc(@(prs) SCALARIZE(@(pars) ...
                                                        WJMODEL_BBFUN(pars, ...
                                                          1, expdat, ...
                                                          K, M, X0, ...
                                                          Fv*Prestress, ...
                                                          L, Txyn, ...
                                                          Qxyn, CFUN, ...
                                                          Npatches, ...
                                                          Nqp, opt), ...
                                                        [prs; ...
                            THETAS_SW(k)], opts), pars0, fopts);  % unconstrained opt
        ERRS_SW(k,:) = WJMODEL_BBFUN(PARS_SW(k,:), 1, expdat, K, M, ...
                                     X0, Fv*Prestress, L, Txyn, Qxyn, ...
                                     CFUN, Npatches, Nqp, opt);
        
        fprintf('Done %d/%d\n', k, length(THETAS_SW));
    end

    save(sprintf('./DATS/ASB_C%s_4IWAN_PSSTIFF_SW_RES.mat', conf), 'PARS_SW', 'err', ...
         'PARMIN_W', 'PARMIN_Z');

    % figure(10)
    % % clf()
    % plot(ERRS_SW(:,1), ERRS_SW(:,2), '.'); hold on
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
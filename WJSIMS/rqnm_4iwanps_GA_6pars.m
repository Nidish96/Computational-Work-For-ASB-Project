clc
clear all 
addpath('../ROUTINES')
addpath('../ROUTINES/WJSIMS')
addpath('../ROUTINES/WJSIMS/IWAN')
addpath('../ROUTINES/FEM')
addpath('../ROUTINES/IWAN')

conf = '1';
Confs = {'1', '2', '3', '4', '4a', 'D1'};

for ci=6:length(Confs)
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

    %% Contact Model
    pars = [9; 20; -0.8; -4; 11.6773; 20];  % [Fs Kt Chi Bt Kp Kn]
    lpars = pars; lpars(copt.lspci) = 10.^(lpars(copt.lspci)); 

    %% Fully Stuck Initial Guess For Prestress Analysis
    Kst = Txyn*diag(reshape([copt.x.T{2}*lpars copt.y.T{2}*lpars ...
                        copt.n.T{1}*lpars]', Npatches*3, 1))*Qxyn;
    K0 = K + Kst;
    X0 = K0\(Prestress*Fv);

    %% Load Experimental Data
    exx = load(sprintf('../EXP_DATA/ASB_%s_20Nm_BB.mat', conf), ...
               'Q', 'W', 'Z');
    Nx = 20; iNs = fix(linspace(1, length(exx.Q), Nx));
    expdat.Q = exx.Q(iNs);
    expdat.W = exx.W(iNs)*2*pi;
    expdat.Z = exx.Z(iNs);
    expdat.D = 2*pi*(expdat.Q.*expdat.W).^2.*expdat.Z;
    
    %% Genetic Algorithm
    gopts = optimoptions('gamultiobj', 'Display', 'iter', ...
                         'PlotFcn', {@gaplotpareto, @gaplotrankhist}, 'MaxGenerations', 100, ...
                         'DistanceMeasureFcn', {@distancecrowding, 'genotype'}, ...
                         'PopulationSize', 100); %,
                                                 %'InitialPopulationMatrix',
                                                 %[PARS]);
    opt = optimoptions('fsolve', 'Display', 'off', ...
                       'SpecifyObjectiveGradient', true); 
    lb = [0, 0, -1, -4, 0, 3];
    ub = [30, 20, 0, 4, 22, 20];
    
    %Genetic with DFUN solution
    [PARS_GA, ERRS_GA] = gamultiobj(@(pars) WJMODEL_BBFUN_MASING(pars, ...
                                                      1, expdat, K, ...
                                                      M, X0, ...
                                                      Fv*Prestress, ...
                                                      L, Txyn, Qxyn, ...
                                                      CFUN, DFUN, ...
                                                      Npatches, opt), ...
                                  length(pars), [], [], [], [], lb, ...
                                  ub, [], gopts);

    save(sprintf('./DATS/ASB_C%s_4IWAN_PSSTIFF_GA_RES.mat', conf), ...
         'PARS_GA', 'ERRS_GA');

    % figure(10)
    % % clf()
    plot(ERRS_GA(:,1), ERRS_GA(:,2), '.-'); hold on
end

%% Backbone function
opt.Display = 'off';

pars = [5.5; 13.8; -0.95; -10; 11.6773; 12.35]
% pars = PARMIN_Z;
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
clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/WJSIMS/')
addpath('../ROUTINES/FEM/')

Confs = {'C1', 'C2', 'C3', 'C4', 'C4a', 'CD1'};

%% Main Configuration
mconf = Confs{6}; % Main configuration
mfname = sprintf('../MATRIX_PREPARE/MAT_NULLRED_%s.mat', mconf);
load(mfname, 'M', 'K', 'L', 'R', 'Fv', 'Th', 'LamT', 'cnum', 'MESH', 'PatchAreas', 'Npatches', 'dofred', 'GTG')

% Experimental Data
exx = load(sprintf('../EXP_DATA/ASB_%s_20Nm_BB.mat', mconf), 'W', 'Z');

% Access Matrices and Functions
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

%% Simulations & Plots
fs = 14;
figure(1);
set(gcf, 'Color', 'w')
clf()
aa = gobjects(size(Confs));

ERRS = cell(size(Confs));
for ic = 1:length(Confs)
    load(sprintf('./DATS/ASB_%s_LINQEP_GA_RES.mat', Confs{ic}), ...
        'PARS_GA', 'ERRS_GA');

    %% Simulations
    ERRS{ic} = zeros(size(ERRS_GA));
    for i=1:size(PARS_GA,1)
        ERRS{ic}(i,:) = QEP_ROOTS(PARS_GA(i,:), K, M*0, M, 1, opts);
    end
    [~, si] = sort(ERRS{ic}(:,1));  ERRS{ic} = ERRS{ic}(si,:);
    

    aa(ic) = loglog(10.^(ERRS{ic}(:,1)/2), 10.^(ERRS{ic}(:,2)/2), '.-', 'LineWidth', 2, 'MarkerSize', 15); hold on
    legend(aa(ic), Confs{ic}, 'fontsize', fs)
    
    if strcmp(Confs{ic}, mconf)
        loglog(10.^(ERRS_GA(:,1)/2), 10.^(ERRS_GA(:,2)/2), 'ko')
    end
end
xlim([1e-16 1e0])

legend(aa(1:end), 'location', 'northwest')

set(gca, 'fontsize', fs);
xlabel('Relative Frequency Error')
ylabel('Relative log-Damping Error')

% export_fig(sprintf('./FIGS/MIXPAR_%s.png', mconf), '-dpng', '-r1200')
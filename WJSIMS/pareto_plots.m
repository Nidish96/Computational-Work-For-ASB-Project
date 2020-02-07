clc
clear all 
addpath('../ROUTINES')
addpath('../ROUTINES/WJSIMS')
addpath('../ROUTINES/WJSIMS/IWAN')
addpath('../ROUTINES/FEM')
addpath('../ROUTINES/IWAN')

Confs = {'C1', 'C2', 'C3', 'C4', 'C4a', 'CD1'};

fs = 14;

%% 4 IWAN Results
figure(1)
set(gcf, 'color', 'w')
clf()

for ci=1:length(Confs)
    conf = Confs{ci};
    load(sprintf('./DATS/ASB_%s_4IWAN_PSSTIFF_GA_RES.mat', conf), ...
        'PARS_GA', 'ERRS_GA');
	ERRS = ERRS_GA;

    [~, si] = sort(ERRS(:,1));
    aa(ci) = loglog(10.^(ERRS(si,1)/2), 10.^(ERRS(si,2)/2), '.-', 'LineWidth', 2); hold on
    legend(aa(ci), conf, 'FontSize', fs)
end
legend(aa(1:end), 'location', 'northwest')

set(gca, 'FontSize', fs)
xlabel('Relative Frequency Error')
ylabel('Relative log-Damping Error')

export_fig('./FIGS/4IWAN_PAR_RES.png', '-dpng', '-r1200')

%% QEP Results
figure(2)
set(gcf, 'color', 'w')
clf()

for ci=1:length(Confs)
    conf = Confs{ci};

    load(sprintf('./DATS/ASB_%s_LINQEP_GA_RES.mat', conf), ...
         'PARS_GA', 'ERRS_GA');
    ERRS = ERRS_GA;
    
    [~, si] = sort(ERRS(:,1));
    aa(ci) = loglog(10.^(ERRS(si,1)/2), 10.^(ERRS(si,2)/2), '.-', 'LineWidth', 2); hold on
    legend(aa(ci), conf, 'FontSize', fs)
end
legend(aa(1:end), 'location', 'northwest')
xlim([1e-16 1e0])

set(gca, 'FontSize', fs)
xlabel('Relative Frequency Error')
ylabel('Relative log-Damping Error')

% export_fig('./FIGS/QEP_PAR_RES.png', '-dpng', '-r1200')
clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/export_fig/')

Confs = {'C1', 'C2', 'C3', 'C4', 'C4a', 'CD1'};
%% Load Data
for i=1:length(Confs)
    Dats(i) = load(sprintf('./ASB_%s_20Nm_BB.mat', Confs{i}));
end

%% Plot Data
fs = 16;

figure(1)
set(gcf, 'Color', 'w')
clf()
aa = gobjects(length(Confs),1);
for i=1:length(Confs)
    aa(i) = semilogx(Dats(i).Q, (Dats(i).W-Dats(i).W(1))/Dats(i).W(1)*100, 'LineWidth', 2); hold on
    legend(aa(i), Confs{i});
end
legend(aa(1:end), 'Location', 'southwest')
xlabel('Modal Amplitude', 'FontSize', fs)
ylabel('Frequency Deviation (%)', 'FontSize', fs)
set(gca, 'FontSize', fs)

export_fig('./FIGS/EXPDAT_W.png', '-dpng', '-r1200')

figure(2)
set(gcf, 'Color', 'w')
clf()
aa = gobjects(length(Confs),1);
for i=1:length(Confs)
    semilogx(Dats(i).Q, (Dats(i).Z-Dats(i).Z(1))/Dats(i).Z(1)*100, 'LineWidth', 2); hold on
end
xlabel('Modal Amplitude', 'FontSize', fs)
ylabel('Damping factor Deviation (%)', 'FontSize', fs)
set(gca, 'FontSize', fs)

export_fig('./FIGS/EXPDAT_Z.png', '-dpng', '-r1200')

disp(table([cellfun(@(c) c(1), {Dats.W})' cellfun(@(c) c(1), {Dats.Z})']))
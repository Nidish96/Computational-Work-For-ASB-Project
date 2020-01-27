clc
clear all 
addpath('../ROUTINES')
addpath('../ROUTINES/WJSIMS')
addpath('../ROUTINES/WJSIMS/IWAN')
addpath('../ROUTINES/FEM')
addpath('../ROUTINES/IWAN')

Confs = {'1', '2', '3', '4', '4a', 'D1'};

figure(1)
clf()

for ci=1:length(Confs)
    conf = Confs{ci};
    load(sprintf('./DATS/ASB_C%s_4IWAN_PSSTIFF_GA_RES.mat', conf), ...
         'PARS_GA', 'ERRS_GA');
    load(sprintf('./DATS/ASB_C%s_4IWAN_PSSTIFF_GA_RES.mat', conf), ...
         'PARS_GA', 'ERRS_GA');    
    
    [~, si] = sort(ERRS_GA(:,1));
    aa(ci) = plot(ERRS_GA(si,1), ERRS_GA(si,2), '.-'); hold on
    legend(aa(ci), ['Configuration ' conf])
end
legend(aa(1:end))
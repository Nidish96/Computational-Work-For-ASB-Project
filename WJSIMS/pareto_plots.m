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
%     try 
%         load(sprintf('./DATS/ASB_C%s_4IWAN_PSSTIFF_SBI_RES.mat', conf), ...
%              'PARS_BI', 'ERRS_BI');    
%         ERRS = ERRS_BI;
%     catch me
%         load(sprintf('./DATS/ASB_C%s_4IWAN_PSSTIFF_GA_RES.mat', conf), ...
%              'PARS_GA', 'ERRS_GA');
%         ERRS = ERRS_GA;
%     end

    load(sprintf('./DATS/ASB_C%s_LINQEP_GA_RES.mat', conf), ...
                'PARS_GA', 'ERRS_GA');
	ERRS = ERRS_GA;
        
    [~, si] = sort(ERRS(:,1));
    aa(ci) = plot(ERRS(si,1)/2, ERRS(si,2)/2, '.-'); hold on
    legend(aa(ci), ['Configuration ' conf])
end
legend(aa(1:end))

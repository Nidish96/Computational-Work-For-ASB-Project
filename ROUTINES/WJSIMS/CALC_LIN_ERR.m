function ERRS_BI = CALC_LIN_ERR(BB, expdat)
% This function recaculates errors given the structure of backbone data and
% experimental data using a linear RMS error metric rather than the
% previous log scale one. Returned frequency errors should be the same. 
%
% Inputs:
% BB - structure of backbones
% expdat - experimental data to compare to
%
% Outputs:
% ERRS_BI - two columns for frequency and damping, rows match BB. 

ERRS_BI = zeros(length(BB), 2);

for jj = 1:length(BB) %loop through all tests
    
    ErrsWZ = [0, 0];
    
    for ii = 1:length(BB{jj}.W)
        
        if isreal(BB{jj}.W(ii))
            ErrsWZ(1) = ErrsWZ(1) + abs((BB{jj}.W(ii)-expdat.W(ii))/expdat.W(ii));
        else
            ErrsWZ(1) = ErrsWZ(1) + 1;
        end
        
        if isreal(BB{jj}.Z(ii))
            ErrsWZ(2) = ErrsWZ(2) + abs((BB{jj}.Z(ii)-expdat.Z(ii))/expdat.Z(ii));
            
            %Previous method
%             ErrsWZ(2) = ErrsWZ(2) + ((log(BB{jj}.Z(ii))-log(expdat.Z(ii)))/log(expdat.Z(ii)))^2;
        else
            ErrsWZ(2) = ErrsWZ(2) + 1;
        end
        
        
    end
    
    %Sqrt the errors to get rms. 
    ERRS_BI(jj, :) = ErrsWZ/length(BB{jj}.Q);
    
end
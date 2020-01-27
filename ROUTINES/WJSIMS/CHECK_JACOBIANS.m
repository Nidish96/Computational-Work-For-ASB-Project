function norm_error = CHECK_JACOBIANS(CFUN, uxyn, pars, normalize)
%CHECK_JACOBIANS numerically computes jacobians of the function to compare
%with analytical jacobians that the function produces. 
% USAGE:
%   norm_error = CHECK_JACOBIANS(CFUN, uxyn, pars, opt);
% INPUTS:
%   uxyn       : Displacements of contact patch, rows are patches, columns
%                are x,y,z
%   pars       : Parameter set for the contact model [Fs/A; Kt/A; Kn/A]
%   opt        : Option matrices for applying parameters to patches and
%                indices on which parameters are log scale. 
%   normalize  : Normalize to the magnitude of the Jacobian or not (Don't
%                   for DFUN)
% OUTPUTS:
%   norm_error : error for the differences between numerical and analytical
%                Jacobian derivative sets. Rows are different perturbation 
%                amounts, columns are 1 = spatial, 2-4 = parameters (x,y,z). This
%                value is normalized to the norm of the analytical Jacobian

    if(~exist('normalize', 'var'))
        normalize = true;
    end
    
    %perturb parameters by these amounts to check for convergence /
    %Jacobian term accuracy. 
    delta_frac = [0.05 0.01 0.005 0.001 0.0005];
    
    norm_error = zeros(length(delta_frac), 4);
    


    [Fxyn, spatial_J, param_J] = CFUN(uxyn, pars);

    if(size(spatial_J, 3) ~= 1)
        %special format for elastic dry frict, reformat:, still doesnt
        %quite give the correct results with thisprogram because this isnt
        %fully correct
        
        spatial_J = [spatial_J(:, 1, 1), spatial_J(:, 3, 3), spatial_J(:, 3, 3)];
    end

    num_param_J = zeros(size(param_J));
    
    %% Loop through delta_frac:
    
    for ii = 1:length(delta_frac)
        
        % %% Disturb Spatially: 
        tic
        
        Fxyn_high = zeros(size(uxyn));
        Fxyn_low = zeros(size(uxyn));
        
        for kk = 1:3 %loop through all dimensions
            
            uxyn_high = uxyn;
            uxyn_low = uxyn;
            
            uxyn_high(:, kk) = uxyn_high(:, kk)*(1+delta_frac(ii));
            uxyn_low(:, kk) = uxyn_low(:, kk)*(1-delta_frac(ii));
            
            Fxyn_high_tmp = CFUN(uxyn_high, pars);
            Fxyn_low_tmp = CFUN(uxyn_low, pars);
            
            Fxyn_high(:, kk) = Fxyn_high_tmp(:, kk);
            Fxyn_low(:, kk) = Fxyn_low_tmp(:, kk);
            
        end

        
        numer_spatial_J = (Fxyn_high - Fxyn_low) ./ (uxyn*2*delta_frac(ii));
        
        norm_error(ii, 1) = norm(numer_spatial_J - spatial_J) / (norm(spatial_J)*normalize + ~normalize);
        
        % %% Disturb the parameters
        
        for jj = 1:length(pars)
            
            %disp(['On parameter ' num2str(jj)]);
            
            high_pars = pars;
            low_pars = pars;
            
            high_pars(jj) = high_pars(jj) * (1+delta_frac(ii));
            low_pars(jj) = low_pars(jj) * (1-delta_frac(ii));
            
            
            Fxyn_high = CFUN(uxyn, high_pars);
            Fxyn_low = CFUN(uxyn, low_pars);
            
            num_partial = (Fxyn_high - Fxyn_low) / (high_pars(jj) - low_pars(jj));
            
            %save x direction
            num_param_J(:,jj,1) = num_partial(:,1);
            
            %save y direction
            num_param_J(:,jj,2) = num_partial(:,2);
            
            %save z direction
            num_param_J(:,jj,3) = num_partial(:,3);
            
            
        end %end of parameter checking loop
        
        norm_error(ii, 2) = norm(num_param_J(:,:,1) - param_J(:,:,1)) / (norm(param_J(:,:,1))*normalize + ~normalize); 
        norm_error(ii, 3) = norm(num_param_J(:,:,2) - param_J(:,:,2)) / (norm(param_J(:,:,2))*normalize + ~normalize);  
        norm_error(ii, 4) = norm(num_param_J(:,:,3) - param_J(:,:,3)) / (norm(param_J(:,:,3))*normalize + ~normalize);

    end %end of perturbation amount loop
    

end
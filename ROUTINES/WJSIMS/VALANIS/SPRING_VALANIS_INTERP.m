function [Fxyn, varargout] = SPRING_VALANIS_INTERP(uxyn, pars, opt, sol)
%BOUC-WEN friction model for contact interface
% USAGE:
%  [Fxyn, varargout] = SPRING_VALANIS(uxyn, pars, opt);
% INPUTS:
%   uxyn       : Displacements of contact patch, rows are patches, columns
%                are x,y,z
%   pars       : Parameter set for the contact model [E0, lambda, Et, kappa, Kn],
%                kappa is the only parameter not normalized to area and not log scale.   
%   opt        : Option matrices for applying parameters to patches and
%                indices on which parameters are log scale. 
%   sol        : Set of ode45 solutions to interpolate on. 
%
% OUTPUTS:
%   Fxyn - Force vectors, rows are partches, columns are x,y,z
%   varargout{1} - Jacobian w.r.t. spatial coordinates
%   varargout{2} - Jacobian w.r.t. parameters

    %convert log parameters to log scale and take derivatives with respect
    %to the conversion. 
    lpars = pars(:);
    lpars(opt.lspci) = 10.^(lpars(opt.lspci));
    dlparsdpars = ones(size(lpars));
    dlparsdpars(opt.lspci) = lpars(opt.lspci).*log(10);    
    
    %split uxy into sign and magnitude:
    suxy = sign(uxyn(:, 1:2));
    uxy = abs(uxyn(:, 1:2));
    
    
    %% Forces    
    Fxyn = zeros(size(uxyn));
    
    % Valanis - X
    % [E0, lambda, Et, kappa, Kn]
%     E0_x = opt.x.T{1}*lpars;
%     lambda_x = opt.x.T{2}*lpars;
%     Et_x = opt.x.T{3}*lpars;
%     kappa_x = opt.x.T{4}*lpars;
    
    %derivatives of the force vector with respect to the valanis model
    %parameters [E0, lambda, Et, kappa]. Pages are X,Y
    dF_dBWpars = zeros(size(uxyn, 1), length(pars)-1, 2);
    
    %spatial derivatives w.r.t x,y
    dF_dxy = zeros(size(uxyn, 1), 2);
    
    for ii = 1:size(uxyn, 1)
            
        %Pull fun_g from inputs
        fun_g = sol{ii}.extdata.odefun;
        
        %get ode45 results:
        if(uxy(ii, 1) <= sol{ii}.x(end))
            %now evaluate at u
            y = deval(sol{ii}, uxy(ii, 1));

            %to match the format that it came out before:
            y = y';
        else
            warning('For patch %d, Increase maximum expected displacement to: %5.6e', ii, uxy(ii, 1))
            
            [~, y] = ode45(fun_g, [sol{ii}.x(end) uxy(ii, 1)], sol{ii}.y(:, end));
        end
        
        Fxyn(ii, 1) = y(end, 1);
        
        %derivatives with respect to the Valanis Parameters 
        dF_dBWpars(ii, :, 1) = y(end, 2:5)';
        
        %calculate the derivative w.r.t. x or y.
        states_dot = fun_g(uxyn(ii, 1), y(end, :)');
        
        dF_dxy(ii, 1) = states_dot(1);
        
        %Fix sign to match with displacement. 
        Fxyn(ii, 1) = suxy(ii, 1) * Fxyn(ii, 1);
        dF_dBWpars(ii, :, 1) = suxy(ii, 1) * dF_dBWpars(ii, :, 1);
        
    end
    
    % Valanis - Y
    % [E0, lambda, Et, kappa, Kn]
%     E0_y = opt.y.T{1}*lpars;
%     lambda_y = opt.y.T{2}*lpars;
%     Et_y = opt.y.T{3}*lpars;
%     kappa_y = opt.y.T{4}*lpars;
    

    for ii = 1:size(uxyn, 1)
        
        %Pull fun_g from inputs
        fun_g = sol{ii}.extdata.odefun;
        
        %Get ode45 results
        if(uxy(ii, 2) <= sol{ii}.x(end))
            %now evaluate at u
            y = deval(sol{ii}, uxy(ii, 2));

            %to match the format that it came out before:
            y = y';
        else
            warning('For patch %d, Increase maximum expected displacement to: %5.6e', ii, uxy(ii, 2))
            
            [~, y] = ode45(fun_g, [sol{ii}.x(end) uxy(ii, 2)], sol{ii}.y(:, end));
        end
        
        Fxyn(ii, 2) = y(end, 1);
        
        %derivatives with respect to the parameters
        dF_dBWpars(ii, :, 2) = y(end, 2:5)';
        
        %calculate the derivative w.r.t. x or y.
        states_dot = fun_g(uxyn(ii, 2), y(end, :)');
        
        dF_dxy(ii, 2) = states_dot(1);
        
        %Fix sign to match with displacement. 
        Fxyn(ii, 2) = suxy(ii, 2) * Fxyn(ii, 2);
        dF_dBWpars(ii, :, 2) = suxy(ii, 2) * dF_dBWpars(ii, :, 2);
    end
    
    
    % Linear Spring - N
    Kn = opt.n.T{1}*lpars;
    Fxyn(:, 3) = Kn.*uxyn(:,3);
    
    
    %if(n_x < 1 | n_y < 1)
    %    disp('n is too small to have real derivatives w.r.t parameters');
    %    disp(['n_x is: ' num2str(n_x(1)) ' n_y is ' num2str(n_y(1))]);
    %end
    
    %% Jacobians
    if nargout>=2
        
        % calculate the set of derivatives: [partial F_x/partial ux, 
        % partial Fy/partial uy, partial Fz/partial uz]
        varargout{1} = [dF_dxy, Kn];
    end
    
    if nargout>=3
        varargout{2} = zeros(size(Fxyn,1), length(lpars), 3);
        
        % Pages are X, Y, Z.
        %Columns are [E0, lambda, Et, kappa, Kn]
        %Rows match displacements still
        
        
        %Valanis - X
        varargout{2}(:,:,1) = [opt.x.T{1}(:, 1).*dF_dBWpars(:, 1, 1), ... % E0
            opt.x.T{2}(:, 2).*dF_dBWpars(:, 2, 1), ... %lambda
            opt.x.T{3}(:, 3).*dF_dBWpars(:, 3, 1), ... %Et
            opt.x.T{4}(:, 4).*dF_dBWpars(:, 4, 1), ... %kappa
            zeros(size(Fxyn,1), 1)]; % Kn
        
        varargout{2}(:,:,1) = varargout{2}(:,:,1).*dlparsdpars';         
        
        %Valanis - Y
        varargout{2}(:,:,2) = [opt.y.T{1}(:, 1).*dF_dBWpars(:, 1, 2), ... % E0
            opt.y.T{2}(:, 2).*dF_dBWpars(:, 2, 2), ... %lambda
            opt.y.T{3}(:, 3).*dF_dBWpars(:, 3, 2), ... %Et
            opt.y.T{4}(:, 4).*dF_dBWpars(:, 4, 2), ... %kappa
            zeros(size(Fxyn,1), 1)]; % Kn
        
        varargout{2}(:,:,2) = varargout{2}(:,:,2).*dlparsdpars';         
        
        
        % Spring - Z
        varargout{2}(:, :, 3) = [zeros(size(uxyn, 1), length(pars)-1) (uxyn(:,3).*opt.n.T{1}(:, 5))].*dlparsdpars';
        
    end
    
    
    
    
    
    %% Define function for use in ode45
    function fun_g = make_g(E0, lambda, Et, kappa)
        %create the function to integrate for ODE45
        %states are: [F (force); partialF/partialE0; partialF/partialLambda;
        %               partialF/partial Et; partial F/partial kappa];
        %It is assumed that sgn(u_dot) = sgn(u) = sgn(F)
        
        fun_g = @valanis_slope;
        
        function states_dot = valanis_slope(u, states)
                   
            %to make all the derivatives readable:
            F = states(1);
            
            %partial of the g function (=dF/du) w.r.t F
            pGpF = (sign(F).*kappa.*lambda.*(((Et.*abs(u)-abs(F)).*lambda)./E0+1))...
                /(((Et.*abs(u)-abs(F)).*kappa.*lambda)./E0+1).^2 ...
                -(sign(F).*lambda)./(((Et.*abs(u)-abs(F)).*kappa.*lambda)./E0+1);
            
            pGpE0 = (((Et*abs(u)-abs(F))*lambda)/E0+1)/(((Et*abs(u)-abs(F))*kappa*lambda)/E0+1)...
                -((Et*abs(u)-abs(F))*lambda)/(E0*(((Et*abs(u)-abs(F))*kappa*lambda)/E0+1))...
                +((Et*abs(u)-abs(F))*kappa*lambda*(((Et*abs(u)-abs(F))*lambda)/E0+1))...
                /(E0*(((Et*abs(u)-abs(F))*kappa*lambda)/E0+1)^2);
            
            pGpLambda = (Et*abs(u)-abs(F))/(((Et*abs(u)-abs(F))*kappa*lambda)/E0+1)...
                -((Et*abs(u)-abs(F))*kappa*(((Et*abs(u)-abs(F))*lambda)/E0+1))...
                /(((Et*abs(u)-abs(F))*kappa*lambda)/E0+1)^2;
            
            pGpEt = (abs(u)*lambda)/(((Et*abs(u)-abs(F))*kappa*lambda)/E0+1)...
                -(abs(u)*kappa*lambda*(((Et*abs(u)-abs(F))*lambda)/E0+1))...
                /(((Et*abs(u)-abs(F))*kappa*lambda)/E0+1)^2;
            
            pGpKappa = -((Et*abs(u)-abs(F))*lambda*(((Et*abs(u)-abs(F))*lambda)/E0+1))...
                /(((Et*abs(u)-abs(F))*kappa*lambda)/E0+1)^2;
            
            states_dot = [E0 .* (1 + lambda./E0.*(Et.*abs(u) - abs(states(1)))) ./ (1 + kappa.*lambda./E0.*(Et.*abs(u) - abs(states(1)))); %dz/du
                         pGpF .* states(2) + pGpE0; %d/du[partial z/ partial E0]                             
                         pGpF .* states(3) + pGpLambda; %d/du[partial z/ partial lambda]       
                         pGpF .* states(4) + pGpEt; %d/du[partial z/ partial Et]
                         pGpF .* states(5) + pGpKappa]; %d/du[partial z/ partial kappa]
            
            %set NaN states to appropriate values assuming kappa=1 for
            %replacement: NOTE: IT IS POSSIBLE THAT IF kappa>1, NaNs can
            %occur and then these expressions are not correct.
            states_dot(isnan(states_dot)) = 0;
            
            states_dot = states_dot.*(kappa~=1) + (kappa==1).*[E0; 1; 0; 0; 0];
                     
             %This block was needed to deal with NaN derivatives for
             %Bouc-Wen:
%             %replace nan states with zero. The only time this should apply
%             %is for z = 0, the last derivative is NaN, but in a limit
%             %should be 0. - For n<1 this is incorrect since the first term
%             %of each derivative goes to zero, but the partialg/partialA =
%             %1.
%             
%             nan_states_dot = isnan(states_dot);
%             states_dot(nan_states_dot) = 0;
%             
%             %updated limit of states_dot for when n<1 and z = 0.
%             states_dot = states_dot + nan_states_dot.*[0; 1; 0; 0]; %if the NaN states include the partial/partial A state, replace with 1. 
%             
        end
    end
 
end
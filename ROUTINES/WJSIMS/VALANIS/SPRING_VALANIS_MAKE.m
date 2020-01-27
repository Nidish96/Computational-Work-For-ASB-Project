function [CFUN, DFUN] = SPRING_VALANIS_MAKE(max_uxy, pars, opt)
%BOUC-WEN friction model for contact interface, creates functions that
%interpolate on data
% USAGE:
%  [CFUN, DFUN] = SPRING_BOUCWEN_MAKE(max_uxy, pars, opt)
% INPUTS:
%   max_uxy    : Maximum displacement to calculate out to. - vector for
%                each patch
%   pars       : Parameter set for the contact model [A eta n Kn],
%                all except for n are normalized to area. 
%   opt        : Option matrices for applying parameters to patches and
%                indices on which parameters are log scale. 
%   
% OUTPUTS:
%   CFUN - function for calculating forces
%   DFUN - function for calculating dissipation
%
%NOTES:
%   Assumes X and Y pars and opts are identical. 
%   Assumes only the outside patches are different, so 1/5 and 2-4 are two
%   sets of identical patches. 

    %convert log parameters to log scale and take derivatives with respect
    %to the conversion. 
    lpars = pars(:);
    lpars(opt.lspci) = 10.^(lpars(opt.lspci));
%     dlparsdpars = ones(size(lpars)); %Not needed in this program
%     dlparsdpars(opt.lspci) = lpars(opt.lspci).*log(10);  %Not needed in this program  
    
    %% Calculate a results set: 
    
    % Valanis - X
    % [E0, lambda, Et, kappa, Kn]
    E0_x = opt.x.T{1}*lpars;
    lambda_x = opt.x.T{2}*lpars;
    Et_x = opt.x.T{3}*lpars;
    kappa_x = opt.x.T{4}*lpars;
    
    %initialize solutions into a cell structure. 
    sol = cell(1,length(E0_x));
    
    
    %just calculate for one set of parameters and assume that force/all results are just
    %scaled by the patch area (in the interp functions)
    for ii = 1:length(E0_x)
        
        %make ode45 function:
        fun_g = make_g(E0_x(ii), lambda_x(ii), Et_x(ii), kappa_x(ii));
        
        %save data for interpolating
        sol{ii} = ode45(fun_g, [0 max_uxy(ii)], zeros(10, 1));
        
    end
    
    %This assumed symmetry is technically incorrect, but produces less
    %error than this method v. the individual patch method, so may be
    %acceptable. 
%     sol{5} = sol{1};
%     sol{3} = sol{2};
%     sol{4} = sol{2};
    
    
    CFUN = @(uxyn, pars) SPRING_VALANIS_INTERP(uxyn, pars, opt, sol);
    DFUN = @(uxyn, pars) SPRING_VALANIS_DISS_INTERP(uxyn, pars, opt, sol);
    

    
    
    %% Define function for use in ode45
    function fun_g = make_g(E0, lambda, Et, kappa)
        %create the function to integrate for ODE45
        %states are: [F (force); partialF/partialE0; partialF/partialLambda;
        %               partialF/partial Et; partial F/partial kappa];
        %It is assumed that sgn(u_dot) = sgn(u) = sgn(F)
        
        fun_g = @valanis_slope;
        
        function states_dot = valanis_slope(u, states)
            
            %Note AF = forward area. 
            
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
                         pGpF .* states(5) + pGpKappa %d/du[partial z/ partial kappa]
                         states(1); %partial AF/partial u = F
                         states(2); %d/du[partial AF/ partial E0] = partial F/partial E0
                         states(3); %d/du[partial AF/ partial lambda] = partial F/partial lambda
                         states(4); %d/du[partial AF/ partial Et] = partial F/partial Et
                         states(5); %d/du[partial AF/ partial kappa] = partial F/partial kappa
                         ];
                         
                         
            %set NaN states to appropriate values assuming kappa=1 for
            %replacement: NOTE: IT IS POSSIBLE THAT IF kappa>1, NaNs can
            %occur and then these expressions are not correct.
            states_dot(isnan(states_dot)) = 0;
            
            states_dot = states_dot.*(kappa~=1) + (kappa==1).*[E0; 1; 0; 0; 0; states(1:5)];
                     
        end
    end
 
end
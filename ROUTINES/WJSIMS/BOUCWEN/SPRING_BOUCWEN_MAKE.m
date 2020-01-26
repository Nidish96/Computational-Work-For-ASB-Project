function [CFUN, DFUN] = SPRING_BOUCWEN_MAKE(max_uxy, pars, opt)
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
    
    % Bouc-Wen - X
    % [A eat n Kn]
    A_x = opt.x.T{1}*lpars;
    eta_x = opt.x.T{2}*lpars;
    n_x = opt.x.T{3}*lpars;
    
    %initialize solutions into a cell structure. 
    sol = cell(1,length(A_x));
    
    
    %just calculate for one set of parameters and assume that force/all results are just
    %scaled by the patch area (in the interp functions)
    for ii = 1:length(A_x)
        
        %vector for preconditioning, do not use this right now since it was
        %not useful before
        pre_vec = ones(8,1);
        
        fun_g = make_g(A_x(ii), eta_x(ii), n_x(ii), pre_vec);
        
        %save data for interpolating
        sol{ii} = ode45(fun_g, [0 max_uxy(ii)], zeros(8, 1));
        
    end
    
    %This assumed symmetry is technically incorrect, but produces less
    %error than this method v. the individual patch method, so may be
    %acceptable. 
%     sol{5} = sol{1};
%     sol{3} = sol{2};
%     sol{4} = sol{2};
    
    
    CFUN = @(uxyn, pars) SPRING_BOUCWEN_INTERP(uxyn, pars, opt, sol);
    DFUN = @(uxyn, pars) SPRING_BOUCWEN_DISS_INTERP(uxyn, pars, opt, sol);
    
    
    %% Define function for use in ode45
    function fun_g = make_g(A, eta, n, pre_vec)
        %create the function to integrate for ODE45
        %states are: [1. z (force); 
        %             2. partialz/partialA; 
        %             3. partialz/partial eta;
        %             4. partial z/partial n
        %             5. AF (Area under forward curve)
        %             6. partial AF/partialA; 
        %             7. partial AF/partial eta
        %             8. partial AF/partial n
        %             ];
        %u and z are assumed positive and assumed to remain positive - this
        %is not actually true anymore. 
        
        fun_g = @bouc_wen_slope;
        
        function states_dot = bouc_wen_slope(~, states)
            
            states_dot = pre_vec.*[A - eta.*abs(states(1)).^n; %dz/du
                         (-eta.*n.*abs(states(1)).^(n-1).*sign(states(1)) ) .* states(2) + 1; %d/du[partial z/ partial A]                             
                         (-eta.*n.*abs(states(1)).^(n-1).*sign(states(1)) ) .* states(3) - abs(states(1)).^n; %d/du[partial z/ partial eta]
                         (-eta.*n.*abs(states(1)).^(n-1).*sign(states(1)) ) .* states(4) - eta.*log(abs(states(1))).*abs(states(1)).^n %d/du[partial z/ partial n]
                         states(1); %partial A/partial u = z
                         states(2); %d/du[partial AF/ partial A] = partial z/partial A
                         states(3); %d/du[partial AF/ partial eta] = partial z/partial eta
                         states(4); %d/du[partial AF/ partial n] = partial z/partial n
                         ]; 
             
            %replace nan states with zero. The only time this should apply
            %is for z = 0, the last derivative is NaN, but in a limit
            %should be 0. - For n<1 this is incorrect since the first term
            %of each derivative goes to zero, but the partialg/partialA =
            %1.
            
            nan_states_dot = isnan(states_dot);
            states_dot(nan_states_dot) = 0;
            
            %updated limit of states_dot for when n<1 and z = 0.
            states_dot = states_dot + nan_states_dot.*[0; 1; 0; 0; 0; 0; 0; 0]; %if the NaN states include the partial/partial A state, replace with 1. 
            
        end
    end
 
end
function [Dxyn, varargout] = SPRING_BOUCWEN_DISS(uxyn, pars, opt)
%BOUC-WEN friction model for contact interface
% USAGE:
%  [Dxyn, varargout] = SPRING_BOUCWEN_DISS(uxyn, pars, opt);
% INPUTS:
%   uxyn       : Displacements of contact patch, rows are patches, columns
%                are x,y,z
%   pars       : Parameter set for the contact model [A eta n Kn],
%                all except for n are normalized to area. 
%   opt        : Option matrices for applying parameters to patches and
%                indices on which parameters are log scale. 
% OUTPUTS:
%   Dxyn - Dissapation in each direction
%   varargout{1} - Jacobian w.r.t. spatial coordinates
%   varargout{2} - Jacobian w.r.t. parameters

    %convert log parameters to log scale and take derivatives with respect
    %to the conversion. 
    lpars = pars(:);
    lpars(opt.lspci) = 10.^(lpars(opt.lspci));
    dlparsdpars = ones(size(lpars));
    dlparsdpars(opt.lspci) = lpars(opt.lspci).*log(10);    
    
    %% Forces
    
    %just evaluate the force on the positive side since the goal is the
    %area which should be the same going either direction
    suxyn = sign(uxyn);
    uxyn = abs(uxyn);
    
    %Don't both evaluating forces by calling the CFUN model like other DISS
    %models, sice they have to come out of ode45 anyway with all derivatives.
    
    %% Forces/Dissipations
    
    Fxyn = zeros(size(uxyn));
    Dxyn = zeros(size(uxyn));
    Areaxyn = zeros(size(uxyn));
    
    % Bouc-Wen - X
    % [A eat n Kn]
    A_x = opt.x.T{1}*lpars;
    eta_x = opt.x.T{2}*lpars;
    n_x = opt.x.T{3}*lpars;
    
    
    %derivatives of the force vector with respect to the bouc-wen
    %parameters [A, alpha, beta, n]. Pages are X,Y
    dF_dBWpars = zeros(size(uxyn, 1), length(pars)-1, 2);
    dArea_dBWpars = zeros(size(uxyn, 1), length(pars)-1, 2);
    
    %spatial derivatives w.r.t x,y
    dFxyndX = zeros(size(uxyn, 1), 2);
    
    
    %loop through all diplacements in this direction
    for ii = 1:size(uxyn, 1)
        
        %vector for preconditioning, do not use this right now since it was
        %not useful before
        pre_vec = ones(8,1);
        
        fun_g = make_g(A_x(ii), eta_x(ii), n_x(ii), pre_vec);
        
        %now evaluate over the region of [0 u], the output u is only
        %maintained for debugging purposes
        [~, y] = ode45(fun_g, [0 uxyn(ii, 1)], zeros(8, 1));
        
%         disp(['ode45 took ' num2str(size(y, 1)) ' steps']);
%         figure; plot(u, y);legend('z', 'partial A', 'partial eta', 'partial n')
        
        Fxyn(ii, 1) = y(end, 1);
        Areaxyn(ii, 1) = y(end, 5);
        
        %derivatives of force with respect to the Bouc-Wen Parameters - [A eta n]
        dF_dBWpars(ii, :, 1) = y(end, 2:4)';
        dArea_dBWpars(ii, :, 1) = y(end, 6:8)';
        
        
        %calculate the derivative w.r.t. x or y.
        states_dot = fun_g(uxyn(ii, 1), y(end, :)');
        
        dFxyndX(ii, 1) = states_dot(1);
        
    end
    
    
    %Applying Masing, dissipation for a full loop.
    Dxyn(:, 1) = 8.*Areaxyn(:, 1) - 4.*Fxyn(:, 1).*uxyn(:, 1);

    
    % Bouc-Wen - Y
    % [A alpha beta n Kn]
    A_y = opt.y.T{1}*lpars;
    eta_y = opt.y.T{2}*lpars;
    n_y = opt.y.T{3}*lpars;
    

    for ii = 1:size(uxyn, 1)
        
        %vector for preconditioning, do not use this right now since it was
        %not useful before
        pre_vec = ones(8,1);
        
        fun_g = make_g(A_y(ii), eta_y(ii), n_y(ii), pre_vec);
        
        %now evaluate over the region of [0 u], the output u is only
        %maintained for debugging purposes
        [~, y] = ode45(fun_g, [0 uxyn(ii, 2)], zeros(8, 1));
        
%         disp(['ode45 took ' num2str(size(y, 1)) ' steps']);
%         figure; plot(u, y);legend('z', 'partial A', 'partial eta', 'partial n')
        
        Fxyn(ii, 2) = y(end, 1);
        Areaxyn(ii, 2) = y(end, 5);
        
        %derivatives with respect to the Bouc-Wen Parameters - [A eta n]
        dF_dBWpars(ii, :, 2) = y(end, 2:4)';
        dArea_dBWpars(ii, :, 2) = y(end, 6:8)';
        
        
        %calculate the derivative w.r.t. x or y.
        states_dot = fun_g(uxyn(ii, 2), y(end, :)');
        
        dFxyndX(ii, 2) = states_dot(1);
        
    end
    
    
    %Applying Masing, dissipation for a full loop.
    Dxyn(:, 2) = 8.*Areaxyn(:, 2) - 4.*Fxyn(:, 2).*uxyn(:, 2);

    
    % Spring - N
    Dxyn(:, 3) = 0;
    
    
    
    %% Jacobians
    if nargout>=2
        
        %columns are derivatives in XYZ directions, rows are patches
        %This expression is correct regardless of friction model given that
        %the correct model is called above for Fxyn ect.
        varargout{1} = [4*Fxyn(:, 1:2)-4*dFxyndX(:, 1:2).*uxyn(:, 1:2), ...
            zeros(size(uxyn(:,3)))];
        
        %correct sign output
        varargout{1} = varargout{1}.*suxyn;
        
    end
    if nargout>=3
        %Rows are patches, columns are parameters, pages are XYZ
        %pars are [A eta Kn]
        varargout{2} = zeros(size(Dxyn,1), length(lpars), 3);
        
        %calculate derivatives of force w.r.t. parameters:
        dFxyndp = zeros(size(Dxyn,1), length(lpars), 3);
        
        %Bouc-Wen - X
        dFxyndp(:,:,1) = [opt.x.T{1}(:, 1).*dF_dBWpars(:, 1, 1), ... % A
            opt.x.T{2}(:, 2).*dF_dBWpars(:, 2, 1), ... %eta
            opt.x.T{3}(:, 3).*dF_dBWpars(:, 3, 1), ... %n
            zeros(size(Fxyn,1), 1)]; % Kn
        
        dFxyndp(:,:,1) = dFxyndp(:,:,1).*dlparsdpars';         
        
        %Bouc-Wen - Y
        dFxyndp(:,:,2) = [opt.y.T{1}(:, 1).*dF_dBWpars(:, 1, 2), ... % A
            opt.y.T{2}(:, 2).*dF_dBWpars(:, 2, 2), ... %eta
            opt.y.T{3}(:, 3).*dF_dBWpars(:, 3, 2), ... %n
            zeros(size(Fxyn,1), 1)]; % Kn
        
        dFxyndp(:,:,2) = dFxyndp(:,:,2).*dlparsdpars';  
        
        %No dissipation/derivatives in Z
        
        % Bouc-wen gen - X
        
        %derivative of forward area:
        varargout{2}(:, :, 1) = [opt.x.T{1}(:, 1).*dArea_dBWpars(:, 1, 1), ... % A
                                opt.x.T{2}(:, 2).*dArea_dBWpars(:, 2, 1), ... %eta
                                opt.x.T{3}(:, 3).*dArea_dBWpars(:, 3, 1), ... %n
                                zeros(size(Fxyn,1), 1)]; % Kn
        
        varargout{2}(:, :, 1) = varargout{2}(:, :, 1).*dlparsdpars';
        
        %convert forward area to dissipation derivative
        varargout{2}(:, :, 1) = -4*dFxyndp(:,:,1).*uxyn(:,1)+8*varargout{2}(:,:,1);

        
        % Bouc-wen gen - Y
        
        %foward area derivative:
        varargout{2}(:, :, 2) = [opt.y.T{1}(:, 1).*dArea_dBWpars(:, 1, 2), ... % A
                                opt.y.T{2}(:, 2).*dArea_dBWpars(:, 2, 2), ... %eta
                                opt.y.T{3}(:, 3).*dArea_dBWpars(:, 3, 2), ... %n
                                zeros(size(Fxyn,1), 1)]; % Kn
        
        varargout{2}(:, :, 2) = varargout{2}(:, :, 2).*dlparsdpars';
        
        %convert forward area to dissipation derivative
        varargout{2}(:, :, 2) = -4*dFxyndp(:,:,2).*uxyn(:,2)+8*varargout{2}(:,:,2);
        
        % Spring - N: No dissipation
        varargout{2}(:, :, 3) = zeros(size(Dxyn,1), length(lpars));
    end
    
    
    
    
    
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
                         states(1); %partial AF/partial u = z
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
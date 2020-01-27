function [Fxyn, varargout] = SPRING_BOUCWEN(uxyn, pars, opt)
%BOUC-WEN friction model for contact interface
% USAGE:
%  [Fxyn, varargout] = SPRING_BOUCWEN(uxyn, pars, opt);
% INPUTS:
%   uxyn       : Displacements of contact patch, rows are patches, columns
%                are x,y,z
%   pars       : Parameter set for the contact model [A eta n Kn],
%                all except for n are normalized to area. 
%   opt        : Option matrices for applying parameters to patches and
%                indices on which parameters are log scale. 
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
    
    
    %% Forces    
    Fxyn = zeros(size(uxyn));
    
    % Bouc-Wen - X
    % [A eat n Kn]
    A_x = opt.x.T{1}*lpars;
    eta_x = opt.x.T{2}*lpars;
    n_x = opt.x.T{3}*lpars;
    
    %derivatives of the force vector with respect to the bouc-wen
    %parameters [A, alpha, beta, n]. Pages are X,Y
    dF_dBWpars = zeros(size(uxyn, 1), length(pars)-1, 2);
    
    %spatial derivatives w.r.t x,y
    dF_dxy = zeros(size(uxyn, 1), 2);
    
    for ii = 1:size(uxyn, 1)
        
        %vector for preconditioning, do not use this right now since it was
        %not useful before
        pre_vec = ones(4,1);
        
        fun_g = make_g(A_x(ii), eta_x(ii), n_x(ii), pre_vec);
        
        %now evaluate over the region of [0 u], the output u is only
        %maintained for debugging purposes
        [~, y] = ode45(fun_g, [0 uxyn(ii, 1)], zeros(4, 1));
        
%         disp(['ode45 took ' num2str(size(y, 1)) ' steps']);
%         figure; plot(u, y);legend('z', 'partial A', 'partial eta', 'partial n')
        
        Fxyn(ii, 1) = y(end, 1);
        
        %derivatives with respect to the Bouc-Wen Parameters - [A eta n]
        dF_dBWpars(ii, :, 1) = y(end, 2:end)';
        
        %calculate the derivative w.r.t. x or y.
        states_dot = fun_g(uxyn(ii, 1), y(end, :)');
        
        dF_dxy(ii, 1) = states_dot(1);
        
    end
    
    % Bouc-Wen - Y
    % [A alpha beta n Kn]
    A_y = opt.y.T{1}*lpars;
    eta_y = opt.y.T{2}*lpars;
    n_y = opt.y.T{3}*lpars;
    

    for ii = 1:size(uxyn, 1)
        
        %vector for preconditioning, do not use this right now since it was
        %not useful before
        pre_vec = ones(4,1);
        
        fun_g = make_g(A_y(ii), eta_y(ii), n_y(ii), pre_vec);
        
        %now evaluate over the region of [0 u], the output u is only
        %maintained for debugging purposes
        [~, y] = ode45(fun_g, [0 uxyn(ii, 2)], zeros(4, 1));
        
%         disp(['ode45 took ' num2str(size(y, 1)) ' steps']);
%         figure; plot(u, y);legend('z', 'partial A', 'partial eta', 'partial n')
        
        Fxyn(ii, 2) = y(end, 1);
        
        %derivatives with respect to the Bouc-Wen Parameters - [A eta n]
        dF_dBWpars(ii, :, 2) = y(end, 2:end)';
        
        %calculate the derivative w.r.t. x or y.
        states_dot = fun_g(uxyn(ii, 2), y(end, :)');
        
        dF_dxy(ii, 2) = states_dot(1);
        
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
        %Columns are [A alpha beta n Kn]
        %Rows match displacements still
        
        
        %Bouc-Wen - X
        varargout{2}(:,:,1) = [opt.x.T{1}(:, 1).*dF_dBWpars(:, 1, 1), ... % A
            opt.x.T{2}(:, 2).*dF_dBWpars(:, 2, 1), ... %eta
            opt.x.T{3}(:, 3).*dF_dBWpars(:, 3, 1), ... %n
            zeros(size(Fxyn,1), 1)]; % Kn
        
        varargout{2}(:,:,1) = varargout{2}(:,:,1).*dlparsdpars';         
        
        %Bouc-Wen - Y
        varargout{2}(:,:,2) = [opt.y.T{1}(:, 1).*dF_dBWpars(:, 1, 2), ... % A
            opt.y.T{2}(:, 2).*dF_dBWpars(:, 2, 2), ... %eta
            opt.y.T{3}(:, 3).*dF_dBWpars(:, 3, 2), ... %n
            zeros(size(Fxyn,1), 1)]; % Kn
        
        varargout{2}(:,:,2) = varargout{2}(:,:,2).*dlparsdpars';         
        
        
        % Spring - Z
        varargout{2}(:, :, 3) = [zeros(size(uxyn, 1), length(pars)-1) (uxyn(:,3).*opt.n.T{1}(:, 4))].*dlparsdpars';
        
    end
    
    
    
    
    
    %% Define function for use in ode45
    function fun_g = make_g(A, eta, n, pre_vec)
        %create the function to integrate for ODE45
        %states are: [z (force); partialz/partialA; partialz/partial eta;
        %              partial z/partial n];
        %u and z are assumed positive and assumed to remain positive - no
        %longer true
        
        fun_g = @bouc_wen_slope;
        
        function states_dot = bouc_wen_slope(~, states)
            
            states_dot = pre_vec.*[A - eta.*abs(states(1)).^n; %dz/du
                         (-eta.*n.*abs(states(1)).^(n-1).*sign(states(1)) ) .* states(2) + 1; %d/du[partial z/ partial A]                             
                         (-eta.*n.*abs(states(1)).^(n-1).*sign(states(1)) ) .* states(3) - abs(states(1)).^n; %d/du[partial z/ partial eta]
                         (-eta.*n.*abs(states(1)).^(n-1).*sign(states(1)) ) .* states(4) - eta.*log(abs(states(1))).*abs(states(1)).^n]; %d/du[partial z/ partial n]
             
            %replace nan states with zero. The only time this should apply
            %is for z = 0, the last derivative is NaN, but in a limit
            %should be 0. - For n<1 this is incorrect since the first term
            %of each derivative goes to zero, but the partialg/partialA =
            %1.
            
            nan_states_dot = isnan(states_dot);
            states_dot(nan_states_dot) = 0;
            
            %updated limit of states_dot for when n<1 and z = 0.
            states_dot = states_dot + nan_states_dot.*[0; 1; 0; 0]; %if the NaN states include the partial/partial A state, replace with 1. 
            
        end
    end
 
end
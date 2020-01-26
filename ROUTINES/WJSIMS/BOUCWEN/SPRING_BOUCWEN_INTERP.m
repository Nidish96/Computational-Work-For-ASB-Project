function [Fxyn, varargout] = SPRING_BOUCWEN_INTERP(uxyn, pars, opt, sol)
%BOUC-WEN friction model for contact interface, interpolates results from
%SPRING_BOUCWEN_MAKE
% USAGE:
%  [Fxyn, varargout] = SPRING_BOUCWEN(uxyn, pars, opt);
% INPUTS:
%   uxyn       : Displacements of contact patch, rows are patches, columns
%                are x,y,z
%   pars       : Just need this for the length of how many parameters there
%                are. 
%   opt        : needed for the calculation of derivatives w.r.t. pars
%   sol        : Set of ode45 solutions to interpolate on. 
%   ode_funs   : Functions for calculating derivatives- not needed,
%                removed.
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
    
%     if( max(max(uxy(:, 1:2))) > max_uxy )
%         error('Error. Maximum initialized displacement exceeded');
%     end
    
    %% Forces    
    Fxyn = zeros(size(uxyn));
    
    % Bouc-Wen - X
%     % [A eat n Kn]
%     A_x = opt.x.T{1}*lpars;
%     eta_x = opt.x.T{2}*lpars;
%     n_x = opt.x.T{3}*lpars;
    
    %derivatives of the force vector with respect to the bouc-wen
    %parameters [A, alpha, beta, n]. Pages are X,Y
    dF_dBWpars = zeros(size(uxyn, 1), length(pars)-1, 2);
    
    %spatial derivatives w.r.t x,y
    dF_dxy = zeros(size(uxyn, 1), 2);
    
    for ii = 1:size(uxyn, 1)
        
        %Pull fun_g from inputs
%         fun_g = ode_funs{ii};
        fun_g = sol{ii}.extdata.odefun;
        
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
        
        %derivatives with respect to the Bouc-Wen Parameters - [A eta n]
        dF_dBWpars(ii, :, 1) = y(end, 2:4)';
        
        %calculate the derivative w.r.t. x or y.
        states_dot = fun_g(uxy(ii, 1), y(end, :)');
        
        dF_dxy(ii, 1) = states_dot(1);
        
        %Fix sign to match with displacement. 
        Fxyn(ii, 1) = suxy(ii, 1) * Fxyn(ii, 1);
        dF_dBWpars(ii, :, 1) = suxy(ii, 1) * dF_dBWpars(ii, :, 1);
        
    end
    
    % Bouc-Wen - Y
%     % [A alpha beta n Kn]
%     A_y = opt.y.T{1}*lpars;
%     eta_y = opt.y.T{2}*lpars;
%     n_y = opt.y.T{3}*lpars;
    

    for ii = 1:size(uxyn, 1)
        
        %Pull fun_g from inputs
%         fun_g = ode_funs{ii};
        fun_g = sol{ii}.extdata.odefun;
        
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
        
        %derivatives with respect to the Bouc-Wen Parameters - [A eta n]
        dF_dBWpars(ii, :, 2) = y(end, 2:4)';
        
        %calculate the derivative w.r.t. x or y.
        states_dot = fun_g(uxy(ii, 2), y(end, :)');
        
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
    

end
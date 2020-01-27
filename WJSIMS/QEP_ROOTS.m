function [s, dsdpar] = QEP_ROOTS(pars, K, C, M, mi, opts)
%QEP_ROOTS returns the roots of and derivatives of the quadratic
%eigenvalue problem.
% USAGE:
%  [s, dsdpar] = QEP_ROOTS(pars, K, M, opts);
%
% INPUTS:
%  pars		:
%  K, M 	:
%  opts 	:
% USAGE:
%  s    	:
%  dsdpar	:
    
    lpars = pars;
    lpars(opts.lspci) = 10.^lpars(opts.lspci);
    dlparsdpars = eye(length(pars));
    dlparsdpars(opts.lspci, opts.lspci) = diag(lpars(opts.lspci))*log(10);
    
    K0 = K;  % Stiffness Matrix
    C0 = C;  % Damping Matrix
    for i=1:length(pars)
        if ~isempty(opts.dK{i})
            K0 = K0 + opts.dK{i}*pars(i);
        end
        if ~isempty(opts.dC{i})
            C0 = C0 + opts.dC{i}*pars(i);
        end
    end
    
    %% Solve QEP using PolyEig
    [V, Z] = polyeig(K0, C0, M);
    [~, si] = sort(abs(Z));
    Z = Z(si);  V = V(:, si);
    
    % Mode of Interest
    s = Z(mi);
    % Vm = V(:, mi);
    Vm = V(:, mi)/sqrt(V(:, mi)'*(M+C0/(2*s))*V(:, mi));  % Mode
                                                          % shape
                                                          % normalization

    dsdpar = zeros(1, length(pars));
    for i=1:length(pars)
        if ~isempty(opts.dK{i})
            % dsdpar(i) = dsdpar(i) + s*(Vm'*opts.dK{i}*Vm)/(Vm'*(K0- ...
            %                                                   s^2*M)*Vm);
            dsdpar(i) = dsdpar(i) + Vm'*opts.dK{i}*Vm;
        end
        if ~isempty(opts.dC{i})
            % dsdpar(i) = dsdpar(i) - s*(Vm'*(s*opts.dC{i})*Vm)/(Vm'* ...
            %                                                   (K0-s^2*M)*Vm);
            dsdpar(i) = dsdpar(i) + Vm'*(s*opts.dC{i})*Vm;
        end        
    end
    
    % %% Solve QEP using Linear recast
    % A = [zeros(size(M)) eye(size(M));
    %      -M\K0 -M\C0];
    % [V, s] = eigs(A, 10, 'SM'); s = diag(s);
    % [~, si] = sort(abs(s));
    % s = s(si);
    % V = V(:, si);
    % V = V./sqrt(diag(V'*eye(size(A))*V)');
    
    % s = s(mi);
    % Vm = V(:, mi);
    
    % dsdpar = zeros(1, length(pars));
    % for i=1:length(pars)
    %     if ~isempty(opts.dK{i})
    %         dsdpar(i) = dsdpar(i) + Vm'*[zeros(size(M,1), size(M,2)*2); ...
    %                             -M\(opts.dK{i}) zeros(size(M))]*Vm;
    %     end
    %     if ~isempty(opts.dC{i})
    %         dsdpar(i) = dsdpar(i) + Vm'*[zeros(size(M,1), size(M,2)*2); ...
    %                             zeros(size(M)) -M\(opts.dC{i})]*Vm;
    %     end        
    % end
    
    dsdpar = dsdpar*dlparsdpars;
end
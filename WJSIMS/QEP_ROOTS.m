function [s, dsdpar] = QEP_ROOTS(pars, K, C, M, mi, varargin)
%QEP_ROOTS returns the roots of and derivatives of the quadratic
%eigenvalue problem.
% USAGE:
%  [s, dsdpar] = QEP_ROOTS(pars, K, M, opts);
%
% INPUTS:
%  pars		: parameters
%  K, C, M 	: Stiffness, Damping, and Mass matrices
%  opts 	: Structure with options
% 	dK, dC		: Stiffness & Damping Jacobian matrices
%	lspci		: Parameter IDs to be considered in
%				log10-scale
%	lspco		: Output IDs to be considered in
%				log10-scale (applicable only when
%				qtty is 'sq_dev' or 'abs_dev')
% 	mode 		: Choose mode of output. One of,
%			  'standard': a+bi
%			  'polar'   : [r; theta]
%			  'sods'    : [Wn; zeta]
%			  (Second-Order-Dynamical-System)
% 	qtty		: Specify quantity to output. One of,
%			  'value'   : value of the root
%			  'abs_dev' : absolute deviation from
%			  		expdat (relative to expdat)
%			  'sq_dev'  : squared deviation from expdat
%			  		(relative to expdat)
% 	expdat		: 2-vector with [W; Z] such that, these can
% 				be compared with the predictions.
% USAGE:
%  s    	: Chosen quantity (root either as a complex number
%  			a+ib, or a 2-vector; or deviation)
%  dsdpar	: Derivatives (either a complex 1x4 vector or a 2x4 matrix)
    opts = struct();
    if nargin==6
        nflds = fieldnames(varargin{1});
        for i=1:length(nflds)
            opts.(nflds{i}) = varargin{1}.(nflds{i});
        end
        if ~ismember('dK', nflds)
            opts.dK = cell(length(pars),1);
        end
        if ~ismember('dC', nflds)
            opts.dC = cell(length(pars),1);
        end
        if ~ismember('lspci', nflds)
            opts.lspci = [];
        end
        if ~ismember('mode', nflds)
            opts.mode = 'standard';
        end
    end
    
    lpars = pars;
    lpars(opts.lspci) = 10.^lpars(opts.lspci);
    dlparsdpars = eye(length(pars));
    dlparsdpars(opts.lspci, opts.lspci) = diag(lpars(opts.lspci))* ...
        log(10);
    
    K0 = K;  % Stiffness Matrix
    C0 = C;  % Damping Matrix
    for i=1:length(lpars)
        if ~isempty(opts.dK{i})
            K0 = K0 + opts.dK{i}*lpars(i);
        end
        if ~isempty(opts.dC{i})
            C0 = C0 + opts.dC{i}*lpars(i);
        end
    end
    
    % %% Solve QEP using PolyEig - Symmetric Case
    % [V, Z] = polyeig(K0, C0, M);
    % [~, si] = sort(abs(Z));
    % Z = Z(si);  V = V(:, si);
    
    % % Mode of Interest
    % s = Z(mi);
    % Vm = V(:, mi);
    % % Amplitude Normalization
    % Vm = V(:, mi)/sqrt(V(:, mi)'*(M+C0/(2*s))*V(:, mi));
    % % Phase Normalization
    % [~,nj] = max(abs(Vm));
    % Vm = Vm*exp(-1j*angle(Vm(nj)));

    % dsdpar = zeros(1, length(lpars));
    % for i=1:length(lpars)
    %     if ~isempty(opts.dK{i})
    %         dsdpar(i) = dsdpar(i) - (Vm'*opts.dK{i}*Vm)/(Vm'*(2*s*M+C0)*Vm);
    %     end
    %     if ~isempty(opts.dC{i})
    %         dsdpar(i) = dsdpar(i) - s*(Vm'*opts.dC{i}*Vm)/(Vm'*(2*s*M+C0)*Vm);
    %     end        
    % end
    
    %% Solve QEP using PolyEig - Asymmetric Case
    [U, Z] = polyeig(K0, C0, M);
    [V, Zp] = polyeig(K0', C0', M');
    
    % Sort
    [~, si] = sort(abs(Z));  Z  = Z(si);   U = U(:, si);
    [~, si] = sort(abs(Zp)); Zp = Zp(si);  V = V(:, si);

    if imag(Z(mi))<=0
        tmp = Z(mi+1);
        Z(mi) = Z(mi+1);
        Z(mi+1) = tmp;
    end
    if imag(Zp(mi))<=0
        tmp = Zp(mi+1);
        Zp(mi) = Zp(mi+1);
        Zp(mi+1) = tmp;        
    end
    % Choose Mode of Interest
%     if abs(abs(Z(mi))-abs(Zp(mi)))>1e-6
%         error('oh no')
%     end
    s = Z(mi);
    Um = U(:, mi);
    Vm = conj(V(:, mi)');

    % % Amplitude Normalization
    % Um = Um./sqrt(Vm*(M+C0/(2*s))*Um);
    % Vm = Vm./sqrt(Vm*(M+C0/(2*s))*Um);
    % % Phase Normalization
    % [~, nj] = max(abs(Um).*abs(Vm'));
    % Um = Um*exp(-1j*angle(Um(nj)));
    % Vm = Vm*exp(-1j*angle(Vm(nj)));
    
    dsdpar = zeros(1, length(lpars));
    for i=1:length(lpars)
        if ~isempty(opts.dK{i})
            dsdpar(i) = dsdpar(i) - (Vm*opts.dK{i}*Um)/(Vm*(2*s*M+C0)*Um);
        end
        if ~isempty(opts.dC{i})
            dsdpar(i) = dsdpar(i) - s*(Vm*opts.dC{i}*Um)/(Vm*(2*s*M+C0)*Um);
        end        
    end
    
    %% Log-Scaling
    dsdpar = dsdpar*dlparsdpars;
    
    %% Convert to appropriate format
    switch (opts.mode)
      case {'standard', 'Standard'}
        s = s;
        dsdpar = dsdpar;
      case {'polar', 'Polar'}
        s = [abs(s); angle(s)];
        dsdpar = [cos(s(2))*real(dsdpar) + sin(s(2))*imag(dsdpar);
                  (cos(s(2))*imag(dsdpar) - sin(s(2))* ...
                   real(dsdpar))/s(1)];
      case {'sods', 'SODS'}
        s = [abs(s); angle(s)];
        dsdpar = [cos(s(2))*real(dsdpar) + sin(s(2))*imag(dsdpar);
                  sin(s(2))*(cos(s(2))*imag(dsdpar) - sin(s(2))* ...
                              real(dsdpar))/s(1)];
        s(2) = -cos(s(2));
    end
    
    if isfield(opts, 'expdat')
        e = (s(:)-opts.expdat(:))./opts.expdat(:);
        dedpar = dsdpar./opts.expdat(:);        
        switch(opts.qtty)
          case 'value'
          case 'abs_dev'
            dsdpar = sign(e).*dedpar;
            s = abs(e);
          case 'sq_dev'
            dsdpar = 2*e.*dedpar;
            s = e.^2;
        end
    end
    
    
    dsdpar(opts.lspco, :) = dsdpar(opts.lspco, :)./s(opts.lspco, :)/log(10);
    s(opts.lspco) = log10(s(opts.lspco));
end
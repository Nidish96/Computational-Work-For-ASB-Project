function [varargout] = MASING_QSMA_CALC(pars, modelmats, expdat, Prestress, opts)
%MASING_QSMA_CALC calculates the hysteresis backbone and estimates the
%modal characteristics given that Masing's assumptions may be made.
%
% USAGE:
%   [Errors] = MASING_QSMA_CALC(pars, modelmats, expdat, Prestress, opts);
%       (or)
%   [Response, Errors] = MASING_QSMA_CALC(pars, modelmats, expdat, Prestress, opts);
% INPUTS:
%   pars            : (Npatches) x (Npars)
%   modelmats       : Structure
%                       M     : Mass matrix
%                       K     : Stiffness matrix
%                       L     : Null-space matrix
%                       Np    : Number of patches
%                       Ftx   : x-tangential force function handle 
%                               [Fx, dFxdux] = Ftx(ux, pars)
%                       Fty   : y-tangential force function handle 
%                               [Fy, dFyduy] = Fty(uy, pars)
%                       Fn    : Normal force function handle 
%                               [Fn, dFndun] = Fn(un, pars)
%   expdat          : Structure
%   Prestress       : Scalar prestress value
%   opts            : Miscellaneous options for solution step
%                       No      : Number of quadrature points for hysteresis
%                                 integration
%                       mdid    : List of eigenmode numbers
%                       md      : Mode to choose from mdid list
%                       minlA   : Minimum Alpha (in logscale)
%                       maxlA   : Maximum Alpha (in logscale)
%                       Na      : Number of Alphas (in logscale)
%                       Display : 2 possible values tested
%                                'iter': display iteration information for
%                                           prestress and progress along
%                                           backbone
%                                'off' : no display
%                       LscPars  : List of parameter columns to be
%                                   interpreted as logscale (10^%T)
%                       SymL     : List for extended symmetry (see below
%                                   for example
%                       PConv    : Function handle for parameter basis
%                                   conversion. Empty if none necessary
%                                       
% OUTPUTS: 2 varargout modes:
%   nargout = 1; Use this for the errors (in log-scale) for Genetic algorithm
%       Errors = log10([Werr; Zerr])
%   nargout = 2; Use this for general response evaluation
%       Outputs = [Q W Z D];
%       Errors = [Werr; Zerr; Derr]
    
    %% Log-Scale Parameters
    Nparsptchs = size(pars, 1);
    if Nparsptchs<modelmats.Np
        tmp = pars;
        pars= zeros(modelmats.Np, size(tmp,2));
        pars(1:Nparsptchs, :) = tmp;
        pars((Nparsptchs+1):end, :) = tmp(opts.SymL, :);
    end
    pars(:, opts.LscPars) = 10.^(pars(:, opts.LscPars));
    if ~isempty(opts.PConv)
        pars = opts.PConv(pars);
    end

    %% Prestress Analysis
    opt = optimoptions('fsolve', 'Display', opts.Display, 'SpecifyObjectiveGradient', true);
    X0 = modelmats.K\(-Prestress*modelmats.Fv);
    
    [Xstat, ~, ~, ~, dRx] = fsolve(@(X) WJRESFUN([X; 0.0], modelmats.K, ...
        X0*0, modelmats.Fv*Prestress, modelmats.L, @(u) modelmats.Ftx(u,pars), ...
        @(u) modelmats.Fty(u,pars), @(u) modelmats.Fn(u,pars), modelmats.Np), X0, opt);
    %% Eigen Analysis
    [V, Dw] = eigs(dRx, modelmats.M, 10, 'SM');
    [Dw, si] = sort(sqrt(abs(diag(Dw)))/(2*pi));
    V = V(:, opts.mdids(opts.md));
    V = V/sqrt(V'*modelmats.M*V);

    %% QSMA
    [La, Wa] = LGWT(opts.No, 0, 1);
    Alphas = logspace(opts.minlA, opts.maxlA, opts.Na);
    La = [La; 1.0];
    
    Q = zeros(opts.Na, 1);
    W = Q;
    Z = Q;
    D = Q;
    Dx = D;
    Dy = D;
    Dn = D;
    Qj = zeros(size(La));
    
    X0 = Xstat + dRx\(Alphas(1)*La(1)*modelmats.M*V);
    opt = optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', true);
    for i=2:opts.Na
        mshape = V;
        for j=1:length(La)
            Xs = fsolve(@(X) WJRESFUN([X; La(j)*Alphas(i)], modelmats.K,...
                modelmats.M*mshape, modelmats.Fv*Prestress, modelmats.L,... 
                @(u) modelmats.Ftx(u,pars), @(u) modelmats.Fty(u,pars),...
                @(u) modelmats.Fn(u,pars), modelmats.Np), X0, opt);
            X0 = Xs;
            Qj(j) = mshape'*modelmats.M*(Xs-Xstat);
            if i>1  % Update mode shape
                mshape = (Xs-Xstat)/(Qj(j));
%                 mshape = mshape;
            end
        end
        Q(i) = Qj(end);
        W(i) = sqrt(Alphas(i)/Q(i));

        % Modal dissipation - inaccurate
%         Harea = Alphas(i)*Q(i) - Alphas(i)*Wa'*Qj(1:end-1);
%         D(i) = -4*Alphas(i)*Q(i) + 8*Harea;
        
%         % Total Dissipation
%         Dx(i) = sum(modelmats.Dtx(modelmats.L(1:3:(modelmats.Np*3),:)*Xs,pars));
%         Dy(i) = sum(modelmats.Dty(modelmats.L(2:3:(modelmats.Np*3),:)*Xs,pars));
%         Dn(i) = sum(modelmats.Dn(modelmats.L(3:3:(modelmats.Np*3),:)*Xs,pars));
        
        % Perturbed Dissipation
        Dx(i) = sum(modelmats.Dtx(modelmats.L(1:3:(modelmats.Np*3),:)*(Xs-Xstat),pars));
        Dy(i) = sum(modelmats.Dty(modelmats.L(2:3:(modelmats.Np*3),:)*(Xs-Xstat),pars));
        Dn(i) = sum(modelmats.Dn(modelmats.L(3:3:(modelmats.Np*3),:)*(Xs-Xstat),pars));
        
        D(i) = Dx(i) + Dy(i) + Dn(i) + 2*pi*(Q(i)*W(i))^2*1e-7*0;
        
        Z(i) = D(i)/(2*pi*(Q(i)*W(i))^2);
        W(i) = W(i)/(2*pi);
        if ~strcmp(opts.Display, 'off')
            fprintf('Done %d/%d\n', i, opts.Na);
        end
    end
    
    %% Errors
    interped = interp1(Q, [W Z D], expdat.Q);
    Werr = rms(interped(:, 1)-expdat.W);
    Zerr = rms((interped(:, 2)-interped(1,2))-(expdat.Z-expdat.Z(1)));
    Derr = rms(interped(:, 3)-expdat.D);
    if sum(isnan([Werr Zerr Derr]))~=0
        disp('Isnan');
%         Werr = rms(expdat.W);
%         Zerr = rms(expdat.Z);
%         Derr = rms(expdat.D);
    end
    
    %% Outputs
    if nargout==1  %% Primary operation for Genetic Algorithm
        varargout{1} = log10([Werr;Zerr]);
    elseif nargout==2  %% Complete Simulation output
        varargout{1} = [Q W Z D];
        varargout{2} = [Werr; Zerr; Derr];
    end
end
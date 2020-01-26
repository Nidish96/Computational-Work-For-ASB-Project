function [opars] = IWANPARCONVERT(ipars, dir)
%IWANPARCONVERT Convert 4-Parameter Iwan Parameters between equivalent
%representations
% USAGE:
%   opars = IWANPARCONVERT(ipars, dir);
% INPUTS:
%   ipars   : Nx4 vector of Iwan Parameters
%   dir     : If 1, ipars=[Fs Kt Chi Beta]
%             If 2, ipars=[R S Chi Phimax]
% OUTPUTS:
%   opars   : Nx4 vector of Iwan Parameters
    opars = ipars;
    switch dir
        case 1
            Fs = ipars(:, 1);  Kt = ipars(:, 2);  Chi = ipars(:, 3); Bt = ipars(:, 4);
            
            Phimax = Fs.*(1+Bt)./(Kt.*(Bt+(Chi+1)./(Chi+2)));
            R = Fs.*(Chi+1)./(Phimax.^(Chi+2).*(Bt+(Chi+1)./(Chi+2)));
            S = Bt.*R.*Phimax.^(Chi+1)./(Chi+1);
            
            opars(:, 1) = R;
            opars(:, 2) = S;
            opars(:, 3) = Chi;
            opars(:, 4) = Phimax;
        case 2
            R = ipars(:, 1);  S = ipars(:, 2);  Chi = ipars(:, 3);  Phimax = ipars(:, 4);
            
            Fs = (S + R.*Phimax.^(Chi+1)./(Chi+2)).*Phimax;
            Kt = S + R.*Phimax.^(Chi+1)./(Chi+1);
            Bt = S./(R.*Phimax.^(Chi+1)./(Chi+1));
            
            opars(:, 1) = Fs;
            opars(:, 2) = Kt;
            opars(:, 3) = Chi;
            opars(:, 4) = Bt;
        otherwise
            error('Input is unknown')
    end
end


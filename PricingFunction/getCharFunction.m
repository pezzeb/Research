function [charQQfunc,charQSfunc] = getCharFunction(modtype,param,S0,K,TTM,r,div)
% GETCHARFUNCTION is a function that give back the characteristic function
%                 for a specfied model.
%
% INPUT:
%   modtype    [string]         : name of the model that we use
%   param      [vector or cell] : parameters for the model
%   S0         [double]         : Price of underlying
%   K          [double]         : Strike price for the option
%   TTM        [double]         : Time to Maturity for the option
%   r          [double]         : continuous interest rate (risk free)
%   div        [double]         : continuous dividend yield
% OUTPUT:
%   charQQfunc [func handle]    : is a [function handle] and is the char-
%                                 acteristic function for the log(S) Random
%                                 variable under Q (bank account as
%                                 numeraire.)
%   charQSfunc [func handle]    : is a [function handle] and is the char-
%                                 acteristic function for the log(S) Random
%                                 variable under Q^S (underlying stock as
%                                 numeraire.)
% *************************************************************************
% List of possible models
% 1.   'BlackScholes'
% 2.   'BSM'                    %USE mu NOT r-q
% 3.   'Heston'
% 4.   'HestonNewParam'
% 5.   'Bates'
% 6.   'BatesNewParam'
% 7.   'BarndorffNielsenShepard'
% 8.   'VG'
% 9.   'CMGY'
% 10.  'CMGYe'
% 11.  'VG-CIR'
% 12.  'VG-GOU'
% 13.  'NIG-CIR'
% 14.  'NIG-GOU'
% 15.  'Merton'
% 16.  'Kou'
% 17.  'CC_'
%    a 'CC_NIG-SA-N
%    b 'CC_NIG-SG-N
%    c 'CC_NIG-IG-N
%    d 'CC_NIG-SIG-N
%    e 'CC_VG-SA-N
%    f 'CC-VG-SG-N
%    g 'CC_VG-IG-N
%    h 'CC_VG-SIG-N
%    i 'CC_CMGY-SA-N
%    j 'CC-CMGY-SG-N
%    k 'CC_CMGY-IG-N
%    l 'CC_CMGY-SIG-N
%    m 'CC_NIG-SA-M
%    n 'CC_NIG-SG-M
%    o 'CC_NIG-IG-M
%    p 'CC_NIG-SIG-M
%    q 'CC_VG-SA-M
%    r 'CC-VG-SG-M
%    s 'CC_VG-IG-M
%    t 'CC_VG-SIG-M
%    u 'CC_CMGY-SA-M
%    v 'CC-CMGY-SG-M
%    x 'CC_CMGY-IG-M
%    y 'CC_CMGY-SIG-M
%    z 'CC_VG-CSA-N


if(length(modtype)>=3 && strcmpi(modtype(1:3),'CC_')) %Carr combinations
    [pX,pZ,pS] = ModelSpecDivide(modtype(4:end));
    modtype    = modtype(1:2);
end

switch lower(modtype)
    case 'blackscholes'
        if(isnumeric(param))
            sigma_bs       = param(1,1);
        elseif(iscell(param))
            sigma_bs       = param{1,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        mu             = r-div;
        x0             = log(S0);
        charQQfunc     = @(phi) exp(1i*(x0 + (mu - 1/2*sigma_bs.^2).*TTM).*phi).*exp(-1/2.*sigma_bs.^2.*phi.^2.*TTM);
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case 'bsm'
        if(isnumeric(param))
            mu       = param(1,1);
            sigma_bs = param(2,1);
        elseif(iscell(param))
            mu       = param{1,1};
            sigma_bs = param{2,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0       = log(S0);
        
        charQQfunc         = @(phi) exp(1i*(x0 + (mu - 1/2*sigma_bs.^2).*TTM).*phi).*exp(-1/2.*sigma_bs.^2.*phi.^2.*TTM);
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case 'heston'
        if(isnumeric(param))
            nu0    = param(1,1);
            kappa  = param(2,1);
            eta    = param(3,1);
            theta  = param(4,1);
            rho    = param(5,1);
        elseif(iscell(param))
            nu0    = param{1,1};
            kappa  = param{2,1};
            eta    = param{3,1};
            theta  = param{4,1};
            rho    = param{5,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0 = log(S0);
        
        d = @(phi) sqrt((rho.*theta.*phi.*1i - kappa).^2 - theta.^2.*(-1i.*phi - phi.^2));
        g = @(phi) (kappa - rho.*theta.*phi.*1i - d(phi))./(kappa - rho.*theta.*phi.*1i + d(phi));
        
        A = @(phi) 1i.*phi.*(x0 + (r - div).*TTM);
        B = @(phi) eta.*kappa./(theta.^2).*( (kappa - rho.*theta.*phi.*1i - d(phi)).*TTM - 2.*log((1 - g(phi).*exp(-d(phi).*TTM))./(1 - g(phi))) );
        C = @(phi) nu0./(theta.^2).*(         kappa - rho.*theta.*1i.*phi - d(phi)).*((1 - exp(-d(phi).*TTM))./(1 - g(phi).*exp(-d(phi).*TTM)));
        
        charQQfunc = @(phi) exp(A(phi) + B(phi) + C(phi));
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case 'hestonnewparam'
        if(isnumeric(param))
            nu0    = param(1,1);
            kappa  = param(2,1);
            etabar = param(3,1);
            Vtheta = param(4,1);
            rho    = param(5,1);
        elseif(iscell(param))
            nu0    = param{1,1};
            kappa  = param{2,1};
            etabar = param{3,1};
            Vtheta = param{4,1};
            rho    = param{5,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0 = log(S0);
        
        d = @(phi) sqrt((rho.*sqrt(Vtheta).*phi.*1i - kappa).^2 - Vtheta.*(-1i.*phi - phi.^2));
        g = @(phi) (kappa - rho.*sqrt(Vtheta).*phi.*1i - d(phi))./(kappa - rho.*sqrt(Vtheta).*phi.*1i + d(phi));
        
        A = @(phi) 1i.*phi.*(x0 + (r - div).*TTM);
        B = @(phi) etabar./(Vtheta).*( (kappa - rho.*sqrt(Vtheta).*phi.*1i - d(phi)).*TTM - 2.*log((1 - g(phi).*exp(-d(phi).*TTM))./(1 - g(phi))) );
        C = @(phi) nu0./(Vtheta).*(         kappa - rho.*sqrt(Vtheta).*1i.*phi - d(phi)).*((1 - exp(-d(phi).*TTM))./(1 - g(phi).*exp(-d(phi).*TTM)));
        
        charQQfunc         = @(phi) exp(A(phi) + B(phi) + C(phi));
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case 'hestoncuiparam'
        warning('MATLAB:s function "integral" can not handle this for some reason')
        if(isnumeric(param))
            nu0    = param(1,1);
            kappa  = param(2,1);
            eta    = param(3,1);
            theta  = param(4,1);
            rho    = param(5,1);
        elseif(iscell(param))
            nu0    = param{1,1};
            kappa  = param{2,1};
            eta    = param{3,1};
            theta  = param{4,1};
            rho    = param{5,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0 = log(S0);
        
        xi = @(phi) kappa - theta.*rho.*1i.*phi;
        d  = @(phi) sqrt((rho.*theta.*phi.*1i - kappa).^2 - theta.^2.*(-1i.*phi - phi.^2)); %this containts xi
        A1 = @(phi) (phi.^2 + 1i.*phi).*sinh(d(phi).*TTM/2);
        A2 = @(phi) d(phi)./nu0.*cosh(d(phi)*TTM/2) + xi(phi)./nu0.*sinh(d(phi)*TTM/2);
        A  = @(phi) A1(phi)./A2(phi);
        D  = @(phi) log(d(phi)./nu0) + (kappa - d(phi))/2 - log((d(phi)+xi(phi))./(2*nu0) + (d(phi)-xi(phi))./(2.*nu0).*exp(-d(phi).*TTM));
        
        charQQfunc = @(phi) exp(1i.*phi.*(x0 + (r - div).*TTM) - kappa.*eta.*rho.*TTM./theta.*1i.*phi - A(phi) + 2*kappa.*eta./(theta.^2).*D(phi));
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case 'bates'
        if(isnumeric(param))
            nu0      = param(1,1);
            kappa    = param(2,1);
            eta      = param(3,1);
            theta    = param(4,1);
            rho      = param(5,1);
            lambda_j = param(6,1);
            mu_j     = param(7,1);
            sigma_j  = param(8,1);
        elseif(iscell(param))
            nu0      = param{1,1};
            kappa    = param{2,1};
            eta      = param{3,1};
            theta    = param{4,1};
            rho      = param{5,1};
            lambda_j = param{6,1};
            mu_j     = param{7,1};
            sigma_j  = param{8,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0 = log(S0);
        
        d = @(phi) sqrt((rho.*theta.*phi.*1i - kappa).^2 - theta.^2.*(-1i.*phi - phi.^2));
        g = @(phi) (kappa - rho.*theta.*phi.*1i - d(phi))./(kappa - rho.*theta.*phi.*1i + d(phi));
        
        A = @(phi) 1i.*phi.*(x0 + (r - div).*TTM);
        B = @(phi) eta.*kappa./(theta.^2).*( (kappa - rho.*theta.*phi.*1i - d(phi)).*TTM - 2.*log((1 - g(phi).*exp(-d(phi).*TTM))./(1 - g(phi))) );
        C = @(phi) nu0./(theta.^2).*(         kappa - rho.*theta.*1i.*phi - d(phi)).*((1 - exp(-d(phi).*TTM))./(1 - g(phi).*exp(-d(phi).*TTM)));
        
        E = @(phi) -lambda_j.*mu_j.*1i.*phi.*TTM   + lambda_j.*TTM.*((1 + mu_j).^(1i*phi).*exp(sigma_j.^2.*(1/2.*1i.*phi).*(1i.*phi - 1)) - 1);
        
        charQQfunc         = @(phi) exp(A(phi) + B(phi) + C(phi) + E(phi));
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case 'batesnewparam'
        if(isnumeric(param))
            nu0      = param(1,1);
            kappa    = param(2,1);
            etabar   = param(3,1);
            Vtheta   = param(4,1);
            rho      = param(5,1);
            lambda_j = param(6,1);
            mu_j     = param(7,1);
            sigma_j  = param(8,1);
        elseif(iscell(param))
            nu0      = param{1,1};
            kappa    = param{2,1};
            etabar   = param{3,1};
            Vtheta   = param{4,1};
            rho      = param{5,1};
            lambda_j = param{6,1};
            mu_j     = param{7,1};
            sigma_j  = param{8,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0 = log(S0);
        
        d = @(phi) sqrt((rho.*sqrt(Vtheta).*phi.*1i - kappa).^2 - Vtheta.*(-1i.*phi - phi.^2));
        g = @(phi) (kappa - rho.*sqrt(Vtheta).*phi.*1i - d(phi))./(kappa - rho.*sqrt(Vtheta).*phi.*1i + d(phi));
        
        A = @(phi) 1i.*phi.*(x0 + (r - div).*TTM);
        B = @(phi) etabar./(Vtheta).*( (kappa - rho.*sqrt(Vtheta).*phi.*1i - d(phi)).*TTM - 2.*log((1 - g(phi).*exp(-d(phi).*TTM))./(1 - g(phi))) );
        C = @(phi) nu0./(Vtheta).*(         kappa - rho.*sqrt(Vtheta).*1i.*phi - d(phi)).*((1 - exp(-d(phi).*TTM))./(1 - g(phi).*exp(-d(phi).*TTM)));
        
        E = @(phi) -lambda_j.*mu_j.*1i.*phi.*TTM   + lambda_j.*TTM.*((1 + mu_j).^(1i*phi).*exp(sigma_j.^2.*(1/2.*1i.*phi).*(1i.*phi - 1)) - 1);
        
        charQQfunc         = @(phi) exp(A(phi) + B(phi) + C(phi) + E(phi));
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case {'barndorffnielsenshepard','bns'}%THIS IS THE GAMMA VERSION - THE IG-OU IS not IMPLEMENTED
        if(isnumeric(param))
            rho    = param(1,1);
            lambda = param(2,1);
            b      = param(3,1);
            a      = param(4,1);
            sigma0 = param(5,1);
        elseif(iscell(param))
            rho    = param{1,1};
            lambda = param{2,1};
            b      = param{3,1};
            a      = param{4,1};
            sigma0 = param{5,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0=log(S0);
        
        f1 = @(phi) 1i.*phi.*rho - 1./lambda.*(phi.^2 + 1i.*phi).*(1 - exp(-lambda.*TTM))./2;
        f2 = @(phi) 1i.*phi.*rho - 1./lambda.*(phi.^2 + 1i.*phi)./2;
        
        A = @(phi) exp(1i.*phi.*(x0 + (r - div - a.*lambda.*rho./(b - rho)).*TTM));
        B = @(phi) exp(-1./lambda.*(phi.^2 + 1i.*phi).*(1 - exp(-lambda.*TTM)).*sigma0.^2./2);
        C = @(phi) exp(a./(b - f2(phi)).*(b.*log((b - f1(phi))./(b - 1i.*phi.*rho)) + f2(phi).*lambda.*TTM));
        
        charQQfunc = @(phi) A(phi).*B(phi).*C(phi);
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case {'errorbns1'}%THIS only implemented for testing purposes, and at the time of writing code I beleive that it is incorrect
        if(isnumeric(param))
            rho    = param(1,1);
            lambda = param(2,1);
            b      = param(3,1);
            a      = param(4,1);
            sigma0 = param(5,1);
        elseif(iscell(param))
            rho    = param{1,1};
            lambda = param{2,1};
            b      = param{3,1};
            a      = param{4,1};
            sigma0 = param{5,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0=log(S0);
        
        f1 = @(phi) 1i.*phi.*rho - 1./1.*(phi.^2 + 1i.*phi).*(1 - exp(-lambda.*TTM))./2;
        f2 = @(phi) 1i.*phi.*rho - 1./1.*(phi.^2 + 1i.*phi)./2;
        
        A = @(phi) exp(1i.*phi.*(x0 + (r - div - a.*lambda.*rho./(b - rho)).*TTM));
        B = @(phi) exp(-1./lambda.*(phi.^2 + 1i.*phi).*(1 - exp(-lambda.*TTM)).*sigma0.^2./2);
        C = @(phi) exp(a./(b - f2(phi)).*(b.*log((b - f1(phi))./(b - 1i.*phi.*rho)) + f2(phi).*lambda.*TTM));
        
        charQQfunc = @(phi) A(phi).*B(phi).*C(phi);
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case {'bns_igou'}
        if(isnumeric(param))
            rho    = param(1,1);
            lambda = param(2,1);
            Gamma  = param(3,1);
            delta  = param(4,1);
            sigma0 = param(5,1);
        elseif(iscell(param))
            rho    = param{1,1};
            lambda = param{2,1};
            Gamma  = param{3,1};
            delta  = param{4,1};
            sigma0 = param{5,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0=log(S0);
        
        f1 = @(phi) 1i.*phi.*rho - 1./lambda.*(phi.^2 + 1i.*phi).*(1 - exp(-lambda.*TTM))./2;
        f2 = @(phi) 1i.*phi.*rho - 1./lambda.*(phi.^2 + 1i.*phi)./2;
        
        A = @(phi) exp(1i.*phi.*(x0 + (r - div - rho.*delta./sqrt(Gamma^2 - 2*rho)).*TTM));
        B = @(phi) exp(-1./lambda.*(phi.^2 + 1i.*phi).*(1 - exp(-lambda.*TTM)).*sigma0.^2./2);
        C = @(phi) exp(delta*(sqrt(Gamma^2 - 2*f1(phi)) - sqrt(Gamma^2 - 2*1i*phi*rho)) + ...
            delta*2*f2(phi)./sqrt(Gamma^2 - 2*f2(phi)).*(atanh(sqrt((Gamma^2 - 2*f1(phi))./(Gamma^2 - 2*f2(phi)))) - atanh(sqrt((Gamma^2 - 2*1i*phi*rho)./(Gamma^2 - 2*f2(phi))))) ...
            );
        
        charQQfunc = @(phi) A(phi).*B(phi).*C(phi);
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case {'errorbns2_1'} %This have all three errors: lambda, tanh/tan, sign
        if(isnumeric(param))
            rho    = param(1,1);
            lambda = param(2,1);
            Gamma  = param(3,1);
            delta  = param(4,1);
            sigma0 = param(5,1);
        elseif(iscell(param))
            rho    = param{1,1};
            lambda = param{2,1};
            Gamma  = param{3,1};
            delta  = param{4,1};
            sigma0 = param{5,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0=log(S0);
        
        f1 = @(phi) 1i.*phi.*rho - 1./1.*(phi.^2 + 1i.*phi).*(1 - exp(-lambda.*TTM))./2;
        f2 = @(phi) 1i.*phi.*rho - 1./1.*(phi.^2 + 1i.*phi)./2;
        
        A = @(phi) exp(1i.*phi.*(x0 + (r - div - rho.*delta./sqrt(Gamma^2 - 2*rho)).*TTM));
        B = @(phi) exp(-1./lambda.*(phi.^2 + 1i.*phi).*(1 - exp(-lambda.*TTM)).*sigma0.^2./2);
        C = @(phi) exp(delta*(sqrt(Gamma^2 - 2*f1(phi)) - sqrt(Gamma^2 - 2*1i*phi*rho)) + ...
            delta*2*f2(phi)./sqrt(2*f2(phi) - Gamma^2).*(atan(sqrt((Gamma^2 - 2*f1(phi))./(2*f2(phi) - Gamma^2))) - atan(sqrt((Gamma^2 - 2*1i*phi*rho)./(2*f2(phi) - Gamma^2)))) ...
            );
        
        charQQfunc = @(phi) A(phi).*B(phi).*C(phi);
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case {'errorbns2_2'} %This have     two   errors: lambda, tanh/tan,     
        if(isnumeric(param))
            rho    = param(1,1);
            lambda = param(2,1);
            Gamma  = param(3,1);
            delta  = param(4,1);
            sigma0 = param(5,1);
        elseif(iscell(param))
            rho    = param{1,1};
            lambda = param{2,1};
            Gamma  = param{3,1};
            delta  = param{4,1};
            sigma0 = param{5,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0=log(S0);
        
        f1 = @(phi) 1i.*phi.*rho - 1./1.*(phi.^2 + 1i.*phi).*(1 - exp(-lambda.*TTM))./2;
        f2 = @(phi) 1i.*phi.*rho - 1./1.*(phi.^2 + 1i.*phi)./2;
        
        A = @(phi) exp(1i.*phi.*(x0 + (r - div - rho.*delta./sqrt(Gamma^2 - 2*rho)).*TTM));
        B = @(phi) exp(-1./lambda.*(phi.^2 + 1i.*phi).*(1 - exp(-lambda.*TTM)).*sigma0.^2./2);
        C = @(phi) exp(delta*(sqrt(Gamma^2 - 2*f1(phi)) - sqrt(Gamma^2 - 2*1i*phi*rho)) + ...
            delta*2*f2(phi)./sqrt(Gamma^2 - 2*f2(phi)).*(atan(sqrt((Gamma^2 - 2*f1(phi))./(Gamma^2 - 2*f2(phi)))) - atan(sqrt((Gamma^2 - 2*1i*phi*rho)./(Gamma^2 - 2*f2(phi))))) ...
            );
        
        charQQfunc = @(phi) A(phi).*B(phi).*C(phi);
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case {'errorbns2_3'} %This have     two   errors: lambda,         , sign
        if(isnumeric(param))
            rho    = param(1,1);
            lambda = param(2,1);
            Gamma  = param(3,1);
            delta  = param(4,1);
            sigma0 = param(5,1);
        elseif(iscell(param))
            rho    = param{1,1};
            lambda = param{2,1};
            Gamma  = param{3,1};
            delta  = param{4,1};
            sigma0 = param{5,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0=log(S0);
        
        f1 = @(phi) 1i.*phi.*rho - 1./1.*(phi.^2 + 1i.*phi).*(1 - exp(-lambda.*TTM))./2;
        f2 = @(phi) 1i.*phi.*rho - 1./1.*(phi.^2 + 1i.*phi)./2;
        
        A = @(phi) exp(1i.*phi.*(x0 + (r - div - rho.*delta./sqrt(Gamma^2 - 2*rho)).*TTM));
        B = @(phi) exp(-1./lambda.*(phi.^2 + 1i.*phi).*(1 - exp(-lambda.*TTM)).*sigma0.^2./2);
        C = @(phi) exp(delta*(sqrt(Gamma^2 - 2*f1(phi)) - sqrt(Gamma^2 - 2*1i*phi*rho)) + ...
            delta*2*f2(phi)./sqrt(2*f2(phi) - Gamma^2).*(atanh(sqrt((Gamma^2 - 2*f1(phi))./(2*f2(phi) - Gamma^2))) - atanh(sqrt((Gamma^2 - 2*1i*phi*rho)./(2*f2(phi) - Gamma^2)))) ...
            );
        
        charQQfunc = @(phi) A(phi).*B(phi).*C(phi);
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case {'errorbns2_4'} %This have     two   errors:       , tanh/tan, sign
        if(isnumeric(param))
            rho    = param(1,1);
            lambda = param(2,1);
            Gamma  = param(3,1);
            delta  = param(4,1);
            sigma0 = param(5,1);
        elseif(iscell(param))
            rho    = param{1,1};
            lambda = param{2,1};
            Gamma  = param{3,1};
            delta  = param{4,1};
            sigma0 = param{5,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0=log(S0);
        
        f1 = @(phi) 1i.*phi.*rho - 1./lambda.*(phi.^2 + 1i.*phi).*(1 - exp(-lambda.*TTM))./2;
        f2 = @(phi) 1i.*phi.*rho - 1./lambda.*(phi.^2 + 1i.*phi)./2;
        
        A = @(phi) exp(1i.*phi.*(x0 + (r - div - rho.*delta./sqrt(Gamma^2 - 2*rho)).*TTM));
        B = @(phi) exp(-1./lambda.*(phi.^2 + 1i.*phi).*(1 - exp(-lambda.*TTM)).*sigma0.^2./2);
        C = @(phi) exp(delta*(sqrt(Gamma^2 - 2*f1(phi)) - sqrt(Gamma^2 - 2*1i*phi*rho)) + ...
            delta*2*f2(phi)./sqrt(2*f2(phi) - Gamma^2).*(atanh(sqrt((Gamma^2 - 2*f1(phi))./(2*f2(phi) - Gamma^2))) - atanh(sqrt((Gamma^2 - 2*1i*phi*rho)./(2*f2(phi) - Gamma^2)))) ...
            );
        
        charQQfunc = @(phi) A(phi).*B(phi).*C(phi);
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);        
    case {'errorbns2_5'} %This have all three errors: lambda,               
        if(isnumeric(param))
            rho    = param(1,1);
            lambda = param(2,1);
            Gamma  = param(3,1);
            delta  = param(4,1);
            sigma0 = param(5,1);
        elseif(iscell(param))
            rho    = param{1,1};
            lambda = param{2,1};
            Gamma  = param{3,1};
            delta  = param{4,1};
            sigma0 = param{5,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0=log(S0);
        
        f1 = @(phi) 1i.*phi.*rho - 1./1.*(phi.^2 + 1i.*phi).*(1 - exp(-lambda.*TTM))./2;
        f2 = @(phi) 1i.*phi.*rho - 1./1.*(phi.^2 + 1i.*phi)./2;
        
        A = @(phi) exp(1i.*phi.*(x0 + (r - div - rho.*delta./sqrt(Gamma^2 - 2*rho)).*TTM));
        B = @(phi) exp(-1./lambda.*(phi.^2 + 1i.*phi).*(1 - exp(-lambda.*TTM)).*sigma0.^2./2);
        C = @(phi) exp(delta*(sqrt(Gamma^2 - 2*f1(phi)) - sqrt(Gamma^2 - 2*1i*phi*rho)) + ...
            delta*2*f2(phi)./sqrt(Gamma^2 - 2*f2(phi)).*(atanh(sqrt((Gamma^2 - 2*f1(phi))./(Gamma^2 - 2*f2(phi)))) - atanh(sqrt((Gamma^2 - 2*1i*phi*rho)./(Gamma^2 - 2*f2(phi))))) ...
            );
        
        charQQfunc = @(phi) A(phi).*B(phi).*C(phi);
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case {'errorbns2_6'} %This have     two   errors:                 , sign
        if(isnumeric(param))
            rho    = param(1,1);
            lambda = param(2,1);
            Gamma  = param(3,1);
            delta  = param(4,1);
            sigma0 = param(5,1);
        elseif(iscell(param))
            rho    = param{1,1};
            lambda = param{2,1};
            Gamma  = param{3,1};
            delta  = param{4,1};
            sigma0 = param{5,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0=log(S0);
        
        f1 = @(phi) 1i.*phi.*rho - 1./lambda.*(phi.^2 + 1i.*phi).*(1 - exp(-lambda.*TTM))./2;
        f2 = @(phi) 1i.*phi.*rho - 1./lambda.*(phi.^2 + 1i.*phi)./2;
        
        A = @(phi) exp(1i.*phi.*(x0 + (r - div - rho.*delta./sqrt(Gamma^2 - 2*rho)).*TTM));
        B = @(phi) exp(-1./lambda.*(phi.^2 + 1i.*phi).*(1 - exp(-lambda.*TTM)).*sigma0.^2./2);
        C = @(phi) exp(delta*(sqrt(Gamma^2 - 2*f1(phi)) - sqrt(Gamma^2 - 2*1i*phi*rho)) + ...
            delta*2*f2(phi)./sqrt(2*f2(phi) - Gamma^2).*(atanh(sqrt((Gamma^2 - 2*f1(phi))./(2*f2(phi) - Gamma^2))) - atanh(sqrt((Gamma^2 - 2*1i*phi*rho)./(2*f2(phi) - Gamma^2)))) ...
            );
        
        charQQfunc = @(phi) A(phi).*B(phi).*C(phi);
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case {'errorbns2_7'} %This have     two   errors:         tanh/tan,     
        if(isnumeric(param))
            rho    = param(1,1);
            lambda = param(2,1);
            Gamma  = param(3,1);
            delta  = param(4,1);
            sigma0 = param(5,1);
        elseif(iscell(param))
            rho    = param{1,1};
            lambda = param{2,1};
            Gamma  = param{3,1};
            delta  = param{4,1};
            sigma0 = param{5,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0=log(S0);
        
        f1 = @(phi) 1i.*phi.*rho - 1./lambda.*(phi.^2 + 1i.*phi).*(1 - exp(-lambda.*TTM))./2;
        f2 = @(phi) 1i.*phi.*rho - 1./lambda.*(phi.^2 + 1i.*phi)./2;
        
        A = @(phi) exp(1i.*phi.*(x0 + (r - div - rho.*delta./sqrt(Gamma^2 - 2*rho)).*TTM));
        B = @(phi) exp(-1./lambda.*(phi.^2 + 1i.*phi).*(1 - exp(-lambda.*TTM)).*sigma0.^2./2);
        C = @(phi) exp(delta*(sqrt(Gamma^2 - 2*f1(phi)) - sqrt(Gamma^2 - 2*1i*phi*rho)) + ...
            delta*2*f2(phi)./sqrt(Gamma^2 - 2*f2(phi)).*(atan(sqrt((Gamma^2 - 2*f1(phi))./(Gamma^2 - 2*f2(phi)))) - atan(sqrt((Gamma^2 - 2*1i*phi*rho)./(Gamma^2 - 2*f2(phi))))) ...
            );
        
        charQQfunc = @(phi) A(phi).*B(phi).*C(phi);
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case 'vg'
        if(isnumeric(param))
            nu    = param(1,1);
            sigma = param(2,1);
            theta = param(3,1);
        elseif(iscell(param))
            nu    = param{1,1};
            sigma = param{2,1};
            theta = param{3,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0 = log(S0);
        
        omega = 1/nu.*log(1 - theta*nu - 1/2*sigma^2*nu);
        
        A = @(phi) exp(1i.*(x0+(r - div + omega).*TTM).*phi);
        B = @(phi) (1 - 1i.*theta.*nu.*phi + 1/2.*sigma.^2.*phi.^2.*nu ).^(-TTM./nu);
        
        charQQfunc = @(phi) A(phi).*B(phi);
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case 'cmgy'
        if(isnumeric(param))
            C      = param(1,1);
            G      = param(2,1);
            M      = param(3,1);
            Y      = param(4,1);
        elseif(iscell(param))
            C      = param{1,1};
            G      = param{2,1};
            M      = param{3,1};
            Y      = param{4,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0=log(S0);
        
        charCGMY   = @(phi) exp(TTM.*C.*gamma(-Y).*( (M - 1i.*phi).^Y - M.^Y + (G+1i.*phi).^Y - G.^Y));
        omega      = -1./TTM.*charCGMY(-1i);% omega      = -1./TTM.*exp(TTM.*C.*gamma(-Y).*( (M-1).^Y - M.^Y + (G+1).^Y - G.^Y));
        
        charQQfunc = @(phi) exp(1i.*phi.*(x0 + (r - div+ omega).*TTM)).*charCGMY(phi);
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case 'cmgye'
        if(isnumeric(param))
            C      = param(1,1);
            G      = param(2,1);
            M      = param(3,1);
            Y      = param(4,1);
            eta    = param(4,1);
        elseif(iscell(param))
            C      = param{1,1};
            G      = param{2,1};
            M      = param{3,1};
            Y      = param{4,1};
            eta    = param{4,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0=log(S0);
        
        charCGMY   = @(phi) exp(TTM.*C.*gamma(-Y).*( (M - 1i.*phi).^Y - M.^Y + (G+1i.*phi).^Y - G.^Y));
        omega      = -1./TTM.*charCGMY(-1i);% omega      = -1./TTM.*exp(TTM.*C.*gamma(-Y).*( (M-1).^Y - M.^Y + (G+1).^Y - G.^Y));
        
        charQQfunc = @(phi) exp(1i.*phi.*(x0 + (r - div + omega - eta.^2/2).*TTM)).*charCGMY(phi).*exp(-phi.^2.*eta.^2/2);
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case 'vg-cir'
        if(isnumeric(param))
            C      = param(1,1);
            G      = param(2,1);
            M      = param(3,1);
            
            kappa  = param(4,1);
            eta    = param(5,1);
            lambda = param(6,1);
            y0     = param(7,1);
        elseif(iscell(param))
            C      = param{1,1};
            G      = param{2,1};
            M      = param{3,1};
            
            kappa  = param{4,1};
            eta    = param{5,1};
            lambda = param{6,1};
            y0     = param{7,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0=log(S0);
        
        charFuncVGL  = @(phi)   ((G.*M)./(G.*M + (M - G).*1i.*phi + phi.^2)).^C;
        %gamma is a function in MATLAB that is used in another switch case and thus to reduce problem this is called gammaF
        gammaF = @(phi) sqrt(kappa.^2 - 2.*lambda.^2.*1i.*phi);
        Nomin1 = exp(kappa.^2.*eta.*TTM./(lambda.^2));
        Nomin2 = @(phi) exp(2.*y0.*1i.*phi./(kappa + gammaF(phi).*coth(gammaF(phi).*TTM/2)));
        Denom1 = @(phi) (cosh(gammaF(phi).*TTM./2) + kappa.*sinh(gammaF(phi).*TTM./2)./gammaF(phi)).^(2.*kappa.*eta./(lambda.^2));
        
        charFuncCIR = @(phi) Nomin1.*Nomin2(phi)./Denom1(phi);
        
        charFuncINN     = charFuncVGL;
        charFuncOUT     = charFuncCIR;
        
        charQQfunc = @(phi) exp(1i.*phi.*((r - div).*TTM + x0)).*(charFuncOUT(-1i.*log(charFuncINN(phi))))./((charFuncOUT(-1i.*log(charFuncINN(-1i)))).^(1i.*phi));
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case 'vg-gou'
        if(isnumeric(param))
            C      = param(1,1);
            G      = param(2,1);
            M      = param(3,1);
            
            lambda = param(4,1);
            a      = param(5,1);
            b      = param(6,1);
            y0     = param(7,1);
        elseif(iscell(param))
            C      = param{1,1};
            G      = param{2,1};
            M      = param{3,1};
            
            lambda = param{4,1};
            a      = param{5,1};
            b      = param{6,1};
            y0     = param{7,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0 = log(S0);
        
        %Time processes
        charFuncVGL  = @(phi)   ((G.*M)./(G.*M + (M - G).*1i.*phi + phi.^2)).^C;
        
        %GAMMA-OU
        Term1 = @(phi) 1i.*phi.*y0./lambda.*(1 - exp(-lambda.*TTM));
        Term2 = @(phi) (lambda.*a)./(1i.*phi - lambda.*b).*(b.*log(b./(b - 1i.*phi./lambda.*(1 - exp(-lambda.*TTM)))) - 1i.*phi.*TTM);
        
        charFuncGOU = @(phi) exp(Term1(phi) + Term2(phi));
        
        charFuncINN     = charFuncVGL;
        charFuncOUT     = charFuncGOU;
        
        charQQfunc = @(phi) exp(1i.*phi.*((r - div).*TTM + x0)).*(charFuncOUT(-1i.*log(charFuncINN(phi))))./((charFuncOUT(-1i.*log(charFuncINN(-1i)))).^(1i.*phi));
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case 'nig-cir'
        if(isnumeric(param))
            alpha  = param(1,1);
            beta   = param(2,1);
            delta  = param(3,1);
            
            kappa  = param(4,1);
            eta    = param(5,1);
            lambda = param(6,1);
            y0     = param(7,1);
        elseif(iscell(param))
            alpha  = param{1,1};
            beta   = param{2,1};
            delta  = param{3,1};
            
            kappa  = param{4,1};
            eta    = param{5,1};
            lambda = param{6,1};
            y0     = param{7,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0 = log(S0);
        
        charFuncNIG = @(phi)    exp(-delta.*(sqrt(alpha.^2 - (beta + 1i.*phi).^2) - sqrt(alpha.^2 - beta.^2)));
        
        gammaF = @(phi) sqrt(kappa.^2 - 2.*lambda.^2.*1i.*phi);
        
        Nomin1 = exp(kappa.^2.*eta.*TTM./(lambda.^2));
        Nomin2 = @(phi) exp(2.*y0.*1i.*phi./(kappa + gammaF(phi).*coth(gammaF(phi).*TTM/2)));
        Denom1 = @(phi) (cosh(gammaF(phi).*TTM./2) + kappa.*sinh(gammaF(phi).*TTM./2)./gammaF(phi)).^(2.*kappa.*eta./(lambda.^2));
        
        charFuncCIR = @(phi) Nomin1.*Nomin2(phi)./Denom1(phi);
        
        charFuncINN     = charFuncNIG;
        charFuncOUT     = charFuncCIR;
        
        charQQfunc = @(phi) exp(1i.*phi.*((r - div).*TTM + x0)).*(charFuncOUT(-1i.*log(charFuncINN(phi))))./((charFuncOUT(-1i.*log(charFuncINN(-1i)))).^(1i.*phi));
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case 'nig-gou'
        if(isnumeric(param))
            alpha  = param(1,1);
            beta   = param(2,1);
            delta  = param(3,1);
            
            lambda = param(4,1);
            a      = param(5,1);
            b      = param(6,1);
            y0     = param(7,1);
        elseif(iscell(param))
            alpha  = param{1,1};
            beta   = param{2,1};
            delta  = param{3,1};
            
            lambda = param{4,1};
            a      = param{5,1};
            b      = param{6,1};
            y0     = param{7,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0 = log(S0);
        
        %Time processes
        charFuncNIG = @(phi)    exp(-delta.*(sqrt(alpha.^2 - (beta + 1i.*phi).^2) - sqrt(alpha.^2 - beta.^2)));
        
        %GAMMA-OU
        Term1 = @(phi) 1i.*phi.*y0./lambda.*(1 - exp(-lambda.*TTM));
        Term2 = @(phi) (lambda.*a)./(1i.*phi - lambda.*b).*(b.*log(b./(b - 1i.*phi./lambda.*(1 - exp(-lambda.*TTM)))) - 1i.*phi.*TTM);
        
        charFuncGOU = @(phi) exp(Term1(phi) + Term2(phi));
        
        charFuncINN     = charFuncNIG;
        charFuncOUT     = charFuncGOU;
        
        charQQfunc = @(phi) exp(1i.*phi.*((r - div).*TTM + x0)).*(charFuncOUT(-1i.*log(charFuncINN(phi))))./((charFuncOUT(-1i.*log(charFuncINN(-1i)))).^(1i.*phi));
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case 'merton'
        if(isnumeric(param))
            lam    = param(1,1);
            mu1    = param(2,1);
            sig    = param(3,1);
            sigs   = param(4,1);
        elseif(iscell(param))
            lam    = param{1,1};
            mu1    = param{2,1};
            sig    = param{3,1};
            sigs   = param{4,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0       = log(S0);
        m        = x0 + (r - div - lam*mu1 - sig^2/2)*TTM;
        mus      = log(1+mu1)-sigs^2/2;
        y        = log(K);
        
        charQQfunc = @(phi) exp(1i.*phi.*m - phi.^2.*sig^2.*TTM/2 - lam.*TTM).*exp(lam.*TTM.*exp(1i.*phi.*mus - phi.^2.*sigs.^2/2));
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case 'kou'
        if(isnumeric(param))
            sig    = param(1,1);
            eta1   = param(2,1);
            eta2   = param(3,1);
            lam    = param(4,1);
            p      = param(5,1);
        elseif(iscell(param))
            sig    = param{1,1};
            eta1   = param{2,1};
            eta2   = param{3,1};
            lam    = param{4,1};
            p      = param{5,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0       = log(S0);
        q        = 1 - p;
        zeta     = p/(eta1-1) - q/(eta2+1);
        
        charQQfunc = @(phi) exp(  1i.*phi.*x0 + TTM.*(  1i.*phi.*(r - div - 1/2*sig.^2 - lam.*zeta) - phi.^2.*sig.^2/2 + lam.*1i.*phi.*(p./(eta1 - 1i.*phi) - q./(eta2 + 1i.*phi) )     )    );
        charQSfunc = @(phi) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+phi);
    case 'cc'
        x0       = log(S0);
        switch upper(pX)
            case 'NIG'
                if(isnumeric(param))
                    nu    = param(1,1);
                    theta = param(2,1);
                    sigma = param(3,1);
                elseif(iscell(param))
                    nu    = param{1,1};
                    theta = param{2,1};
                    sigma = param{3,1};
                else
                    error('Integralpricing:Unknown variable type')
                end
                
                psiX = @(u) sigma*(nu/sigma - sqrt(  nu^2/(sigma^2)-2*theta/(sigma^2)*1i*u + u.^2) );
            case 'VG'
                if(isnumeric(param))
                    G     = param(1,1);
                    M     = param(2,1);
                    C     = param(3,1);
                elseif(iscell(param))
                    G     = param{1,1};
                    M     = param{2,1};
                    C     = param{3,1};
                else
                    error('Integralpricing:Unknown variable type')
                end
                
                psiX = @(u) C*log(   (G*M)./( G*M+(M-G)*1i*u + u.^2   ));
            case 'CGMY'
                if(isnumeric(param))
                    Cp    = param(1,1);
                    M     = param(2,1);
                    G     = param(3,1);
                    Yp    = param(4,1);
                    Ym    = param(5,1);
                    zeta  = param(6,1);                    
                elseif(iscell(param))
                    Cp    = param{1,1};
                    M     = param{2,1};
                    G     = param{3,1};
                    Yp    = param{4,1};
                    Ym    = param{5,1};
                    zeta  = param{6,1};                     
                else
                    error('Integralpricing:Unknown variable type')
                end
                %gamma here is the function gamma and NOT a variable
                psiX = @(u) Cp.*(gamma(-Yp).*((M-1i*u).^Yp - M.^Yp) + zeta.*gamma(-Ym).*( (G+1i*u).^Ym - G.^Ym) );
            otherwise
                error('USERdef: The Model: %s, that is requested is not implemented',modtype)
        end
        switch upper(pZ)
            case 'SA'
                if(isnumeric(param))
                    kappa = param(1,2);
                    eta   = param(2,2);
                    lambda= param(3,2);
                    y0    = param(4,2);
                elseif(iscell(param))
                    kappa = param{1,2};
                    eta   = param{2,2};
                    lambda= param{3,2};
                    y0    = param{4,2};
                end
                gam       = @(u) sqrt(kappa^2 - 2*lambda^2*1i*u); %gamma is a function in MATLAB
                phiZ      = @(u) exp((kappa^2*eta*TTM)/(lambda^2))./((cosh(gam(u)*TTM/2) + kappa./gam(u).*sinh(gam(u).*TTM/2)).^(2*kappa*eta/(lambda^2))).*... %A(t,u)
                    exp(( (2*1i*u)./(kappa + gam(u).*coth(gam(u)*TTM/2)) )*y0);
            case 'SG'
                if(isnumeric(param))
                    kappa = param(1,2);
                    lambda= param(2,2);
                    xi    = param(3,2);
                    rho   = param(4,2);
                    y0    = param(5,2);
                elseif(iscell(param))
                    kappa = param{1,2};
                    lambda= param{2,2};
                    xi    = param{3,2};
                    rho   = param{4,2};
                    y0    = param{5,2};
                end
                psi_U     = @(u) 1i*u*lambda./(1/xi - 1i*u);
                
                INTE      = @(a,b) log( (b + a*(1-exp(-kappa*TTM))/kappa + 1i/xi).^( lambda./(kappa - 1i*xi*(a + kappa*b))) .* (a + kappa*b - kappa*(b + a*(1-exp(-kappa*TTM))/kappa)).^( (lambda*(a + kappa*b)*xi)./( kappa*((a + kappa*b)*xi+1i*kappa) ) ) ) - ...
                                   log( (b                               + 1i/xi).^( lambda./(kappa - 1i*xi*(a + kappa*b))) .* (a + kappa*b - kappa*(b                              )).^( (lambda*(a + kappa*b)*xi)./( kappa*((a + kappa*b)*xi+1i*kappa) ) ) );
                
                Phi_t     = @(a,b) exp(1i*a*y0*(1-exp(-kappa*TTM))/(kappa)).*exp(INTE(a,b));
            case 'IG'
                if(isnumeric(param))
                    kappa = param(1,2);
                    nu    = param(2,2);
                    rho   = param(3,2);
                    y0    = param(4,2);
                elseif(iscell(param))
                    kappa = param{1,2};
                    nu    = param{2,2};
                    rho   = param{3,2};
                    y0    = param{4,2};
                end
                psi_U     = @(u) nu - sqrt(nu^2 - 2*1i*u);
                
                INTE      = @(a,b) (2*(sqrt(nu^2 - 2*1i*(b + a*(1-exp(-kappa*TTM))/kappa) ))/kappa - 2/(kappa^(3/2)) * (nu^2*kappa - 2*1i*(a+kappa*b))./(sqrt(nu^2*kappa - 2*1i*(a + kappa*b))).*atanh( (sqrt(kappa)*sqrt(nu^2-2*1i*(b + a*(1-exp(-kappa*TTM))/kappa)) )./(sqrt( nu^2*kappa - 2*1i*(a + kappa*b))) ) - nu*log(a + kappa*b - kappa*(b + a*(1-exp(-kappa*TTM))/kappa))/kappa) - ...
                                   (2*(sqrt(nu^2 - 2*1i*(b                              ) ))/kappa - 2/(kappa^(3/2)) * (nu^2*kappa - 2*1i*(a+kappa*b))./(sqrt(nu^2*kappa - 2*1i*(a + kappa*b))).*atanh( (sqrt(kappa)*sqrt(nu^2-2*1i*(b                              )) )./(sqrt( nu^2*kappa - 2*1i*(a + kappa*b))) ) - nu*log(a + kappa*b - kappa*(b                              ))/kappa);
                
                Phi_t     = @(a,b) exp(1i*a*y0*(1-exp(-kappa*TTM))/(kappa)).*exp(INTE(a,b));
            case 'SIG'
                if(isnumeric(param))
                    kappa = param(1,2);
                    nu    = param(2,2);
                    rho   = param(3,2);
                    y0    = param(4,2);
                elseif(iscell(param))
                    kappa = param{1,2};
                    nu    = param{2,2};
                    rho   = param{3,2};
                    y0    = param{4,2};
                end
                psi_U = @(u) 1i*u./sqrt(nu^2 - 2*1i*u);
                
                INTE      = @(a,b) ((sqrt(nu^2 - 2*1i*(b + a*(1-exp(-kappa*TTM))/kappa) ))/kappa - (2*1i*(a + kappa*b))./(kappa^(3/2)*sqrt(nu^2*kappa - 2*1i*(a + kappa*b))).*atanh(   (sqrt(kappa)*sqrt(nu^2-2*1i*(b + a*(1-exp(-kappa*TTM))/kappa)) )./(sqrt( nu^2*kappa - 2*1i*(a + kappa*b))) )) - ...
                                   ((sqrt(nu^2 - 2*1i*(b                              ) ))/kappa - (2*1i*(a + kappa*b))./(kappa^(3/2)*sqrt(nu^2*kappa - 2*1i*(a + kappa*b))).*atanh(   (sqrt(kappa)*sqrt(nu^2-2*1i*(b                              )) )./(sqrt( nu^2*kappa - 2*1i*(a + kappa*b))) ));
                
                Phi_t     = @(a,b) exp(1i*a*y0*(1-exp(-kappa*TTM))/(kappa)).*exp(INTE(a,b));
            case 'CSA'
                if(~strcmpi(pX,'VG'))
                    if(~strcmpi(pS,'N'))
                        error('USERdef: the model X that you want to use: %s is not valid with CSA. Only ''VG'' is valid. The architecture you want to use %s is not valid in CSA. Only N is valid.',pX,pS)
                    else
                        error('USERdef: the model X that you want to use: %s is not valid with CSA. Only ''VG'' is valid.',pX)
                    end
                else
                    if(strcmpi(pS,'N'))
                        %CORRECT - do nothing
                    else
                        error('The architecture you want to use %s is not valid in CSA. Only N is valid.',pS)
                    end
                end
                if(isnumeric(param))
                    kappa = param(1,2);
                    eta   = param(2,2);
                    lambda= param(3,2);
                    rho2  = param(4,2);
                    y0    = param(5,2);
                elseif(iscell(param))
                    kappa = param{1,2};
                    eta   = param{2,2};
                    lambda= param{3,2};
                    rho2  = param{4,2};
                    y0    = param{5,2};
                end
                gam       = @(a) sqrt(kappa^2 - 2*lambda^2*1i*a); %gamma is a function in MATLAB
                Phi_t     = @(a,b) exp((kappa^2*eta*TTM)/(lambda^2))./((cosh(gam(a)*TTM/2) + (kappa-1i*b*lambda^2)./gam(a).*sinh(gam(a).*TTM/2)).^(2*kappa*eta/(lambda^2))).*... %A(t,u)
                                   exp(( (1i*b.*(gam(a).*cosh(gam(a)*TTM/2) - kappa*sinh(gam(a)*TTM/2)) + 2*1i*a.*sinh(gam(a)*TTM/2))./( gam(a).*cosh(gam(a)*TTM/2) + (kappa-1i*b*lambda^2).*sinh(gam(a)*TTM/2)  ) )*y0);
            otherwise
                error('USERdef: The Model: %s, that is requested is not implemented',modtype)
        end
        
        %Put the parts together
        switch upper(pZ)
            case 'SA'
                switch upper(pS)
                    case 'N'
                        charQQfunc = @(u) exp(1i*u*(x0+(r-div)*TTM)).*phiZ( -1i*psiX(u) )./( (phiZ( -1i*psiX(-1i) )).^(1i*u) );
                        charQSfunc = @(u) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+u);
                    case 'M'
                        charQQfunc = @(u) exp(1i*u*(x0+(r-div)*TTM)).*phiZ( -1i*psiX(u) - u.*psiX(-1i));
                        charQSfunc = @(u) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+u);
                    otherwise
                        error('USERdef: The Model: %s, that is requested is not implemented',modtype)
                end
            case {'SG','IG','SIG'}
                switch pS
                    case 'N'
                        charQQfunc = @(u) exp(1i*u*(x0+(r-div)*TTM)).*Phi_t(-1i*psiX(u),rho*u).*exp(-1i*u*log(Phi_t(-1i*psiX(-1i),-1i*rho)));
                        charQSfunc = @(u) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+u);
                    case 'M'
                        charQQfunc = @(u) exp(1i*u*(x0+(r-div)*TTM -psi_U(-1i*rho)*TTM)).*Phi_t( -1i*psiX(u) - u.*psiX(-1i),rho*u);
                        charQSfunc = @(u) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+u);
                    otherwise
                        error('USERdef: The Model: %s, that is requested is not implemented',modtype)
                end
            case 'CSA'
                charQQfunc = @(u) exp(1i*u*(x0+(r-div)*TTM)).*Phi_t(-1i*psiX(u),rho2*u).*exp(-1i*u.*log(Phi_t(-1i*psiX(-1i),-1i*rho2 )));
                charQSfunc = @(u) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+u);
            otherwise
                error('USERdef: The Model: %s, that is requested is not implemented',modtype)
        end
        
    case 'fmls'    
        if(isnumeric(param))
            alp    = param(1,1);
            sig    = param(2,1);
        elseif(iscell(param))
            alp    = param{1,1};
            sig    = param{2,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0       = log(S0);
        mu       = sig^alp*sec(pi*alp/2);
        
        charQQfunc = @(u) exp(1i*u*(x0 + (r - div + mu)*TTM) - TTM*(1i*u*sig).^alp*sec(alp/2*pi));
        charQSfunc = @(u) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+u);
    case 'fmls+'    
        if(isnumeric(param))
            alp    = param(1,1);
            sig1   = param(2,1);
            sig2   = param(3,1);
        elseif(iscell(param))
            alp    = param{1,1};
            sig1   = param{2,1};
            sig2   = param{3,1};
        else
            error('Integralpricing:Unknown variable type')
        end
        
        x0       = log(S0);
        
        charQQfunc = @(u) exp(1i*u*(x0 + (r - div + sig1^alp*sec(pi*alp/2) - sig2^2)*TTM) - TTM*(1i*u*sig1).^alp.*sec(pi*alp/2) - TTM*(u*sig2).^2);
        charQSfunc = @(u) exp(-(r-div).*TTM-log(S0)).*charQQfunc(-1i+u);
    otherwise
        error('USERdef: The Model: %s, that is requested is not implemented',modtype)
end

end
%% HELP FUNCTIONS
function [p1,p2,p3] = ModelSpecDivide(modelname)
pos  = strfind(modelname,'-');
p1   = modelname(1:pos(1)-1);
p2   = modelname(pos(1)+1:pos(2)-1);
p3   = modelname(pos(2)+1:end);
end
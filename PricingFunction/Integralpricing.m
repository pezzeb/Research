function [call_price,put_price] = Integralpricing(modtype,inttype,param,S0,K,TTM,r,div,alpha)
%INTEGRALPRICING is a function that uses integrals for pricing Call and Put
%                options.
%
% INPUT:
%   modtype    [string]         : name of the model that we use
%   inttype    [string]         : name of the integral method that we use
%   param      [vector or cell] : parameters for the model
%   S0         [double]         : Price of underlying
%   K          [double]         : Strike price for the option
%   TTM        [double]         : Time to Maturity for the option
%   r          [double]         : continuous interest rate (risk free)
%   div        [double]         : continuous dividend yield
%   alpha      [double]         : the damping constant (used in some methods)
% OUTPUT:
%   call_price [double]         : price of the Call option
%   put_price  [double]         : price of the Put  option

if(8==nargin)
    alpha=1.25;
end

[charQQfunc,charQSfunc] = getCharFunction(modtype,param,S0,K,TTM,r,div);

switch lower(inttype)
    case 'alpha'                     %Carr Madan                                                                    
        ModifiedCharFunc = @(phi) exp(-r*TTM).*charQQfunc(phi - (alpha+1)*1i)./(alpha^2 + alpha - phi.^2 + 1i*(2*alpha +1).*phi);
        IntegrandeN = @(phi) real(exp(-1i.*phi.*log(K)).*ModifiedCharFunc(phi));
        intValue = integral(IntegrandeN,0,inf);
        
        call_price = exp(-log(K)*alpha)/pi*intValue;
        put_price = call_price + K.*exp(-r*TTM)- S0.*exp(-div*TTM);
    case 'not_alpha'                 %Kahl Jackels pricing formula Not-so-complex Logarithms                        
        integrand1 = @(u)              real( exp(-1i.*u.*log(K)).*charQQfunc(u - 1i)./(1i.*u));
        integrand2 = @(u) K.*          real( exp(-1i.*u.*log(K)).*charQQfunc(u)./     (1i.*u));
        
        int1 = integral(integrand1,0,inf);
        int2 = integral(integrand2,0,inf);
        
        call_price = 1/2*(S0*exp(-div*TTM) - exp(-r*TTM)*K) + exp(-r*TTM)/pi*(int1 - int2);
        put_price = call_price + K.*exp(-r*TTM)- S0.*exp(-div*TTM);
    case 'not_alpha_probcomp'        %Probability Computation                                                       
        integrand1 = @(u) real( exp(-1i.*u.*log(K)).*charQSfunc(u)./(1i.*u));
        integrand2 = @(u) real( exp(-1i.*u.*log(K)).*charQQfunc(u)./(1i.*u));
        
        probQS = 1/2 + 1/pi*integral(integrand1,0,inf);
        probQQ = 1/2 + 1/pi*integral(integrand2,0,inf);
        
        call_price = S0.*exp(-div.*TTM).*probQS - K.*exp(-r.*TTM).*probQQ;
        put_price  = call_price + K.*exp(-r*TTM)- S0.*exp(-div*TTM);
    case 'alpha_psquad'              %Carr Madan ps quad implementation                                             
        ModifiedCharFunc = @(phi) exp(-r*TTM).*charQQfunc(phi - (alpha+1)*1i)./(alpha^2 + alpha - phi.^2 + 1i*(2*alpha +1).*phi);
        IntegrandeN = @(phi) real(exp(-1i.*phi.*log(K)).*ModifiedCharFunc(phi));
        intValue = PSquadmethod(IntegrandeN,0,inf,'quadgk');
        
        call_price = exp(-log(K)*alpha)/pi*intValue;
        put_price = call_price + K.*exp(-r*TTM)- S0.*exp(-div*TTM);
    case 'not_alpha_psquad'          %Kahl Jackels pricing formula Not-so-complex Logarithms ps quad implementation 
        integrand1 = @(u)              real( exp(-1i.*u.*log(K)).*charQQfunc(u - 1i)./(1i.*u));
        integrand2 = @(u) K.*          real( exp(-1i.*u.*log(K)).*charQQfunc(u)./     (1i.*u));
        
        int1 = PSquadmethod(integrand1,0,inf,'quadgk');
        int2 = PSquadmethod(integrand2,0,inf,'quadgk');
        
        if(int1==inf || int2==inf)
            call_price=inf;
            put_price = inf;
        else
            call_price = 1/2*(S0*exp(-div*TTM) - exp(-r*TTM)*K) + exp(-r*TTM)/pi*(int1 - int2);
            put_price = call_price + K.*exp(-r*TTM)- S0.*exp(-div*TTM);
        end
    case 'not_alpha_probcomp_psquad' %Probability Computation ps quad implementation                                
        integrand1 = @(u) real( exp(-1i.*u.*log(K)).*charQSfunc(u)./(1i.*u));
        integrand2 = @(u) real( exp(-1i.*u.*log(K)).*charQQfunc(u)./(1i.*u));
        
        probQS = 1/2 + 1/pi*PSquadmethod(integrand1,0,inf,'quadgk');
        probQQ = 1/2 + 1/pi*PSquadmethod(integrand2,0,inf,'quadgk');
        
        call_price = S0.*exp(-div.*TTM).*probQS - K.*exp(-r.*TTM).*probQQ;
        put_price  = call_price + K.*exp(-r*TTM)- S0.*exp(-div*TTM);
end

end
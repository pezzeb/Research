function [Oc_out,Op_out] = price_options(Oc,Op,Qd,K,T,rf,type_price,AD_scalar)
%PRICE_OPTIONS Price call and put options
%   Price options without interpolation

if(strcmp(type_price,'interpolation'))
    [Oc_out,Op_out] = interpol_pricing(Oc,Op,Qd,K,T,rf,AD_scalar);
elseif(strcmp(type_price,'exact_match'))
    
    delta_T = f_Delta_T(T);
    discount = delta_T.*rf(1:end-1);
    
    Oc_out = price_single_option(Oc,Qd,K,T,discount,'call');
    Op_out = price_single_option(Op,Qd,K,T,discount,'put');
    
else
    error('ERROR in option pricing');
end
end

function [Ox_out] = price_single_option(Ox,Qd,S,T,discount,type)

Ox_out = revADm(zeros(size(Ox,1),1));

if(strcmp(type,'call'))
    for i_op=1:size(Ox,1)
        
        %Option information
        str_op = Ox(i_op,1);
        ttm_op = Ox(i_op,2);
        
        %other information
        payoff  = max(S-str_op,0);         % payoff vector
        ttm_id = find(T==ttm_op);   % id for the time to maturity
        qd = Qd(:,ttm_id);          % the distribution for the TTM
        
        Ox_out(i_op) = exp(cum_disc(discount,ttm_id-1))*payoff'*qd;
    end
else
    for i_op=1:size(Ox,1)
        
        %Option information
        str_op = Ox(i_op,1);
        ttm_op = Ox(i_op,2);
        
        %other information
        payoff  = max(str_op-S,0);         % payoff vector
        ttm_id = find(T==ttm_op);   % id for the time to maturity
        qd = Qd(:,ttm_id);          % the distribution for the TTM
        display('minustecken som saknas? i prissättninge på diskfaktorn')
        Ox_out(i_op) = exp(cum_disc(discount,ttm_id-1))*payoff'*qd;
    end
end

end
function out = cum_disc(discount,id)

out = sum(discount(1:id));

end
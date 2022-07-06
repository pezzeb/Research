function [Oc_out,Op_out] = price_options_old(Oc,Op,Qd,S,T,rf)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

delta_T = f_Delta_T(T);
discount = delta_T.*rf(1:end-1);

Oc_out = price_single_option(Oc,Qd,S,T,discount);
Op_out = price_single_option(Op,Qd,S,T,discount);
end
function [Ox_out] = price_single_option(Ox,Qd,S,T,discount)

Ox_out = [Ox,zeros(size(Ox,1),1)];

for i_op=1:size(Ox,1)
    
    %Option information
    str_op = Ox(i_op,1);
    ttm_op = Ox(i_op,2);
    
    %other information
    payoff  = S-str_op;         % payoff vector
    ttm_id = find(T==ttm_op);   % id for the time to maturity
    qd = Qd(:,ttm_id);          % the distribution for the TTM
    
    Ox_out(i_op,4) = cum_disc(discount,ttm_id)*payoff'*qd; 
end

end
function out = cum_disc(discount,id)

out = sum(discount(1:id));

end
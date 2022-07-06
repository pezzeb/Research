function [ FROM, TO ] = f_FROM_TO(type,mu,sigma,T,K,epsilon_min,epsilon_max,Delta_c,fix_up,fix_down)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

FROM = zeros(length(K),length(T)-1);
TO = zeros(size(FROM));



if(strcmp(type,'math'))
Delta_T = f_Delta_T(T);
S_u = max(K);
S_l = min(K);

for k = 1:length(K);
    for t = 1:length(T)-1
        
        mu_t    = mu(k,t)*Delta_T(t);
        sigma_t = sqrt(Delta_T(t))*sigma(k,t);
        
        r_min = mu_t + sigma_t*norminv(epsilon_max,0,1);
        r_max = mu_t + sigma_t*norminv(1-epsilon_min,0,1);

        %The r_max and r_min is fix - and therefore they should not write
        %to the tape
        b_u = ceil(r_max/Delta_c);
        b_l = ceil(-r_min/Delta_c);

        
        FROM(k,t)   = max(k-b_l  ,1);
        TO(k,t)     = min(k + b_u,length(K));
    end
end
elseif(strcmp(type,'fix'))
    K_max = length(K);
    for k = 1:length(K);
        for t = 1:length(T)-1
            FROM(k,t)   = max(k - fix_down,1);
            TO(k,t)     = min(k + fix_up,K_max);
        end
    end
end


end


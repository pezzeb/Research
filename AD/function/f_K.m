function [K_out,i_ini] = f_K(type,Delta_c,Sc,Su,Sl,Ko)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if(strcmp(type,'equdistant_NO_option'))
    n_u = ceil(log(Su/Sc)/Delta_c);
    n_l = ceil(log(Sc/Sl)/Delta_c);
    K_out = Sc*exp((-n_l:1:n_u)*Delta_c)';
    i_ini = find(Sc==K_out);
elseif(strcmp(type,'equidistant'))
    n_u = ceil(log(Su/Sc)/Delta_c);
    n_l = ceil(log(Sc/Sl)/Delta_c);
    Kc = Sc*exp((-n_l:1:n_u)*Delta_c)';
    K_out = sort([Kc;Ko]);
    
elseif(strcmp(type,'math'))
    Ko = sort(Ko);
    
    %Multiple entires
    if(Su<max(Ko) || Sl>min(Ko))
        error('Su and Sl does NOT cover the range op option prices');
    elseif(sum(Ko==Sc)==1)
        Kt=[Sl;Ko;Su];
    else
        [v,i_ini_t] = sort([Ko;Sc]);
        Kt = [Sl;v;Su];
        n_ini = length(v);
        i_new = find(i_ini_t == n_ini);
    end
    
    K_out(1,1) = Sl;
    
    for t=1:length(Kt)-1
       
        n_inc = floor(log(Kt(t+1)/Kt(t))/Delta_c);       
        if (0==n_inc)
            K_out(end+1,1) = Kt(t+1);
        else
            Delta_inc = log(Kt(t+1)/Kt(t))/n_inc;
            
            for i_inc=1:n_inc
                K_out(end+1,1) = K_out(end,1)*exp(Delta_inc);
            end

        end
        if(t ==i_new)
            i_ini = length(K_out);
        end
    end
end
end

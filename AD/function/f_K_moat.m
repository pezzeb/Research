function [K_with_moat,I_moat,i_ini_out] = f_K_moat(type,K,lower_moat,upper_moat,i_ini_in)


if(strcmp(type,'simple'))
    
    if(isempty(lower_moat))
        %THE LOWER MOAT IS EMPTY
        n_l = 0;
    else
        if(max(lower_moat)<min(K))
            %LOWER IS OK
            n_l = length(lower_moat);
            i_ini_out = i_ini_in + n_l;
        else
            error('The LOWER MOAT is not consitent with the grid it should be applied to');
        end
    end
    
    if(isempty(upper_moat))
        %THE UPPER MOAT IS EMPTY
        n_u = 0;
    else
        if(min(upper_moat)>max(K))
            %UPEPR IS OK
            n_u = length(upper_moat);
        else
            error('The UPPER MOAT is not consitent with the grid it should be applied to');
        end
    end

    n   = length(K);
    K_with_moat = sort([upper_moat;K;lower_moat]);
    I_moat = [(1:1:n_l)';((n_l+n)+1:1:(n_l+n+n_u))'];    
else
    error('Not implemented functionallity');
end

end
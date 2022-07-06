function [step_length_out] = step_length(u_bar_val,Delta_u,lowest_allowed_vol,line_search_type,step_len_percentage)

if(strcmp(line_search_type,'fix'))
    step_length_out = 1/100;
elseif(strcmp(line_search_type,'simple'))
    step_length_out = 1;
    while(min(min(u_bar_val+step_length_out*Delta_u))<lowest_allowed_vol)
       step_length_out = step_length_out/2; 
    end
    step_length_out = step_len_percentage*step_length_out;
elseif(strcmp(line_search_type,'optimal'))
    temp_step       = (u_bar_val-lowest_allowed_vol)./Delta_u;
    temp_step_gz    = temp_step(logical((temp_step<0).*(temp_step>-1)));
    
    biggest = -max(temp_step_gz);
    
    step_length_out = step_len_percentage*biggest;
    if(isempty(step_length_out))
        step_length_out = 1;
    end
end

end
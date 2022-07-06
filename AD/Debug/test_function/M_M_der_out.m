function [ out ] = M_M_der_out(u,v,constant_factor)
%UNTITLED4 Summary of this function goes here

[n_row_u,n_col_u] = size(u);
[n_row_v,n_col_v] = size(v);
out = cell(n_row_u,n_col_v);

if (strcmp(constant_factor,'first'))
    n_el_v = numel(v);
    
    for i_r=1:n_row_u
        for i_c=1:n_col_v
            out{i_r,i_c} = first(u,v,i_r,i_c,n_el_v,n_row_u,n_col_u);
        end
    end
elseif (strcmp(constant_factor,'second'))
    
    n_el_u = numel(u);
    
    for i_r=1:n_row_u
        for i_c=1:n_col_v
            out{i_r,i_c} = second(u,v,i_r,i_c,n_el_u,n_row_u,n_col_u);
        end
    end
elseif(strcmp(constant_factor,'none'))
    n_el_u = numel(u);
    n_el_v = numel(v);
    
    for i_r=1:n_row_u
        for i_c=1:n_col_v
            
            %NOTICE THE PART2 correponds to that the second is non constant
            
            part2=  first(u,v,i_r,i_c,n_el_v,n_row_u,n_col_u);
            part1= second(u,v,i_r,i_c,n_el_u,n_row_u,n_col_u);
            out{i_r,i_c} = [part1 part2];
        end
    end
    
end

end

function temp = first(u,v,i_r,i_c,n_el,n_row_u,n_col_u)

v = u(i_r,:);
temp = zeros(1,n_el);
temp((i_c-1)*n_col_u+1:i_c*n_col_u) = v;

end

function temp = second(u,v,i_r,i_c,n_el,n_row_u,n_col_u)

idx = i_r:n_row_u:n_el;
temp = zeros(1,n_el);
temp(idx) = v(:,i_c);

end
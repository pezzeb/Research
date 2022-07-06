function [b_transpose,H,nabla_g ] = reverse_procedure(n_in_var,gc,gp,D,b_e,H_h,b_h_transpose )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


g_id = [gc.id;gp.id];
g_val= [gc.val;gp.val];
nabla_g = zeros(length(g_id),n_in_var);

for i_nabla=1:length(g_id)
    nabla_g(i_nabla,:) = reverse_tape(n_in_var,g_id(i_nabla));
end

H_g = nabla_g'*D*nabla_g;
b_g_transpose = -(b_e-g_val)'*D*nabla_g;

b_transpose = b_g_transpose + b_h_transpose;
H           = H_h + H_g;


end


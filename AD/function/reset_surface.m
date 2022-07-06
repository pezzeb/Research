function [ mu_hat, u_bar, mu_surf, U, u_bar_val, U_hat, U_scale, U_id] = reset_surface(n_in_var,U_val,diff_Tc,r_f,div_yield)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

reset_tape(n_in_var);
U        = f_Pi_ad(U_val);
U_id     = U.id;
U_scale  = repmat(sqrt([diff_Tc;0]'),size(U,1),1);
U_hat    = U.*U_scale;
u_bar    = vec_v(U);
u_bar_val= vec_v(U_val);

%TÄNK PÅ ATT VI PADDAR MED NOLL FÖR ATT FÅ DIMENSIONERNA GÅR IHOP
mu_surf = f_mu(U,r_f,div_yield);
mu_hat    = mu_surf.*repmat([diff_Tc;0]',size(mu_surf,1),1);

end


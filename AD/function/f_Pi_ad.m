function Pi = f_Pi_ad(Pi_val)
global tape
n_in_var = numel(Pi_val);
Pi = revADm(Pi_val,reshape(1:1:n_in_var,size(Pi_val)));
end
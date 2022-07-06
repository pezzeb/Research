function [diff_iv_c, diff_iv_p, diff_pr_c, diff_pr_p] = IV_difference_InS( gc,gp,Oc,Op,IV_market_c,IV_market_p,Sini,r_fcon,div_con)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

gc_val = gc.val;
gp_val = gp.val;

IV_surf_c = zeros(length(gc_val),1);
IV_surf_p = zeros(length(gp_val),1);

%Call options
for i=1:length(gc_val)
    try
        IV_surf_c(i,1) = calcBSImpVol(1 ,gc_val(i),Sini,Oc(i,1),Oc(i,2),r_fcon(1),div_con(1));
    catch me
        IV_surf_c(i,1) = NaN;
    end
end

for i=1:length(gp_val)
    try
        IV_surf_p(i,1) = calcBSImpVol(-1,gp_val(i),Sini,Op(i,1),Op(i,2),r_fcon(1),div_con(1));
    catch me
        IV_surf_p(i,1) = NaN;
    end
end

diff_iv_c = IV_surf_c-IV_market_c;
diff_iv_p = IV_surf_p-IV_market_p;

diff_pr_c = gc_val-Oc(:,3);
diff_pr_p = gp_val-Op(:,3);
end


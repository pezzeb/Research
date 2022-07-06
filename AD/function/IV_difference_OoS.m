function [diff_cell] = IV_difference_OoS( n_out_of_sample_series,cell_IV_market_out,cell_out_of_sample,cell_out_of_sample_price,Sini,r_fcon,div_con)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

diff_cell = cell(n_out_of_sample_series,2);
for i_oos_iv=1:n_out_of_sample_series
    
    IV_market_c_out     = cell_IV_market_out{i_oos_iv,1};
    IV_market_p_out     = cell_IV_market_out{i_oos_iv,2};
    
    Price_c             = cell_out_of_sample_price{i_oos_iv,1};
    Price_p             = cell_out_of_sample_price{i_oos_iv,2};
    
    Oc_out          = cell_out_of_sample{i_oos_iv,1};
    Op_out          = cell_out_of_sample{i_oos_iv,2};
    
    diff_c = [];
    diff_p = [];
    for j=1:size(Oc_out,1)
        try    
            IV_surf_c_out    = calcBSImpVol( 1, Price_c(j), Sini, Oc_out(j,1), Oc_out(j,2), r_fcon(1),div_con(1));
            diff_c(end+1,1)  = IV_surf_c_out - IV_market_c_out(j);
        catch me
            diff_c(end+1,1)  = NaN;
        end
        end
    
    for k=1:size(Op_out,1)
        try
            IV_surf_p_out    = calcBSImpVol(-1, Price_p(k), Sini, Op_out(k,1), Op_out(k,2), r_fcon(1),div_con(1));
            diff_p(end+1,1)  = IV_surf_p_out - IV_market_p_out(k);
        catch me
            diff_p(end+1,1)  = NaN;
        end
    end
    diff_cell{i_oos_iv,1}      = diff_c;
    diff_cell{i_oos_iv,2}      = diff_p;
end
end


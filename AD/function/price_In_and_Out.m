function [ gc,gp,cell_out_of_sample_price ] = price_In_and_Out(Oc,Op,Qd,K,Tc,r_f,type_price,n_out_of_sample_series,cell_out_of_sample)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


%In sample
[gc,gp] = price_options(Oc,Op,Qd,K,Tc,r_f,type_price,'AD');
%Out of sample
if(0<n_out_of_sample_series)
    for i_oos_price=1:n_out_of_sample_series
        [gc_out,gp_out] = price_options(cell_out_of_sample{i_oos_price,1},cell_out_of_sample{i_oos_price,2},Qd,K,Tc,r_f,type_price,'scalar');
        cell_out_of_sample_price{i_oos_price,1} = gc_out;
        cell_out_of_sample_price{i_oos_price,2} = gp_out;
    end
else
    cell_out_of_sample_price = {};
    cell_out_of_sample_price = {};
end
end


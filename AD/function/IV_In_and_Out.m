function [IV_market_c,IV_market_p,cell_IV_market_out] = IV_In_and_Out(Oc,Op,Sini,n_out_of_sample_series,cell_out_of_sample)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

for k =1:size(Oc(:,3),1)
    try
        IV_market_c(k,1) = calcBSImpVol(1 ,Oc(k,3),Sini,Oc(k,1),Oc(k,2),0,0);
    catch me
        IV_market_c(k,1) = NaN;
        display('Problem in market')
    end
end
for k =1:size(Op(:,3),1)
    try
        IV_market_p(k,1) = calcBSImpVol(-1,Op(k,3),Sini,Op(k,1),Op(k,2),0,0);
    catch me
        IV_market_p(k,1) = NaN;
        display('Problem in market')
    end
end


%Out of Sample
cell_IV_market_out = cell(n_out_of_sample_series,2);
for i=1:n_out_of_sample_series
    cell_IV_tempc = [];
    cell_IV_tempp = [];
    
    Oc_out = cell_out_of_sample{i,1};
    Op_out = cell_out_of_sample{i,2};
    
    for j=1:length(Oc_out(:,3))
        try
            cell_IV_tempc(j,1) = calcBSImpVol(1 ,Oc_out(j,3),Sini,Oc_out(j,1),Oc_out(j,2),0,0);
            cell_IV_tempp(j,1) = calcBSImpVol(-1,Op_out(j,3),Sini,Op_out(j,1),Op_out(j,2),0,0);
        catch meiday
        end
    end
    cell_IV_market_out{i,1} = cell_IV_tempc;
    cell_IV_market_out{i,2} = cell_IV_tempp;
end


end


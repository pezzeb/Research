function [ out ] = M_v_der_out(u,v)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[n_row_u,n_col_u] = size(u);
n_el = numel(u);

out = cell(n_row_u,1);

for i=1:n_row_u
    temp = zeros(1,n_el);
    temp(i:n_row_u:n_el) = v';
    out{i,1} = temp;
end

end


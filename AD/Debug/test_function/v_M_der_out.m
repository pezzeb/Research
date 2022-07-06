function [ out ] = v_M_der_out(u,v)
%UNTITLED4 Summary of this function goes here
%   v is a row vector

[n_row_u,n_col_u] = size(u);
n_el = numel(u);

out = cell(1,n_col_u);

for i=1:n_col_u
    temp = zeros(1,n_el);
    temp((i-1)*n_row_u+1:i*n_row_u) = v;
    out{1,i} = temp;
end

end


function [ str_out ] = der_from_matrix(X_val)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[n_row,n_col] = size(X_val);

for i=1:n_row
    for j=1:n_col
        
        T = zeros(n_row,n_col);
        T(i,j) = X_val(i,j);
        str_out{i,j} = T(:)';
        
    end
end

end


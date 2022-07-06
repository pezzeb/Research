function [ rev_out ] = reverse_range( n_in_var,id_range )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
length(id_range)
pause()
for iter=1:length(id_range)
    iter
    rev_out(iter,:) = reverse_tape(n_in_var,id_range(iter));
end


end


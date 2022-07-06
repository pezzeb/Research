function [ r_mat ] = f_r_mat(K)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

r_mat = zeros(length(K),length(K));

for i=1:length(K)
    r_mat(i,:) = log(K/K(i))';
end

end


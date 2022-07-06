function [ mu ] = f_mu(sigma,rf,delta,risk_aversion)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if(3==nargin)
    risk_aversion = 0;
end

[n_row,n_col] = size(sigma);
 
if(size(rf,1) == n_col)
    mu = (1-risk_aversion)*sigma.^2 + repmat(rf',n_row,1) - delta;
else
    %ERROR
end

end
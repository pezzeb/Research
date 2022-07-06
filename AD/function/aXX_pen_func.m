function [ aXX ] = aXX_pen_func(len_K,len_T,pen_func,ref_pen)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

t=1:1:len_T;
k=1:1:len_K;

aXX  = ref_pen * ones(len_K,len_T);
aXX  = aXX    .* pen_func(k,t);

end
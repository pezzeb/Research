function [ test_vector ] = get_tv( n,name,element )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

test_vector = cell(0,0);    %ALLOCATION

if(length(n)==length(name) && length(name)==length(element))
    len_n = length(n);

for i=1:len_n
    
   test_vector(length(test_vector)+1) = name{i};
   
   temp_element = element{i};
   
   test_vector(length(test_vector)+1:length(test_vector)+n(i)) = temp_element(:)';
    
end
else
   error('ERROR i TEST'); 
end


end


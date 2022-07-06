function [start_surf] = start_solution(K,T,sigma,type)


if (strcmp('constant',type))
    start_surf = ones(length(K),length(T))*sigma;
else
    %ERROR
end

end
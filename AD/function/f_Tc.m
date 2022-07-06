function [Tc_out] = f_Tc(type_in,Delta_t,T_max,T_abs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if(strcmp(type_in,'relative_equdistant'))
    Tc_out = (0:Delta_t:T_max)';
elseif(strcmp(type_in,'relative_NONequdistant'))
    Tc_out(1,1) = 0;
    for kk=1:length(Delta_t)
        Tc_out(kk+1,1) = Tc_out(kk,1) + Delta_t(kk);
    end
elseif(strcmp(type_in,'absolute'))
    Tc_out = T_abs;
end
end

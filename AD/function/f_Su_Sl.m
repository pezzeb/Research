function [ Su,Sl ] = f_Su_Sl(epsilon_u,epsilon_l,Sc,mu,sigma)
%F_SU_SL Returns the Su and Sl values.
%   Detailed explanation goes here

Su = Sc*exp(mu+sigma*norminv(1-epsilon_u,0,1));
Sl = Sc*exp(mu+sigma*norminv(epsilon_l,0,1));

end


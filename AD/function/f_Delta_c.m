function [ Delta_c_out ] = f_Delta_c(type,epsilon,mu,sigma,Delta_c)
%F_E Returns a matrix E.
%   PARAMETERS:
%   Delta_c  - is a user defined Delta_c that is choosen as default
%   p_limits - is the probability limits for delta_c

if(strcmp('user',type));
    Delta_c_out = Delta_c;
elseif(strcmp('math',type))
    
    done = false;
    Delta_c_out = 5;
    while (not(done))
        k = ceil(1/2*(ceil(2*mu/Delta_c_out)-1));
        done = (normcdf(Delta_c_out/2+(2*k+1)*Delta_c_out/2,mu,sigma) - ...
                normcdf(Delta_c_out/2+(2*k-1)*Delta_c_out/2,mu,sigma))<epsilon;
        if(done)
            break
        else
           Delta_c_out = Delta_c_out/2; 
        end
    end
    
end

end


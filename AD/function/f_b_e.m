function [rho] = f_b_e(Oc,Op)
%F_RHo Creates rho from to option objects.
%   The function takes two option objects a returns the market prices for
%   the call and put options. The results is an column vector with the call
%   options prices at the top. (TESTED)

if(nargin == 1)
    Op=[];
end
if(isempty(Oc) && isempty(Op))
    rho = [];
elseif(not(isempty(Oc)) && not(isempty(Op)))
    rho = [Oc(:,3);Op(:,3)];
elseif(isempty(Op))
    rho = [Oc(:,3)];
else
    rho = [Op(:,3)];
end
    
end
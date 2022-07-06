function [q,q_exp] = f_node_q(from,to,r_mat_row,r_f,mu_hat_rev,sigma_hat_rev,n_m,sigma_scale)
global q_exceptions_neg;
global q_exceptions_nan;
%F_NODE_Q Returns the distribution for a single node
%   Detailed explanation goes here

mu_hat      = forADm(mu_hat_rev.val,2*sigma_hat_rev.val*sigma_scale);
sigma_hat   = forADm(sigma_hat_rev.val,sigma_scale);          

% mu_hat      = mu_hat_rev;
% sigma_hat   = sigma_hat_rev;

%% Init
if(6==nargin)
    n_m = 2;
end
r_mat_temp = r_mat_row(from:to)';

mc = normal_central_moments(sigma_hat,n_m);

p_bar_not_norm = normpdfSPEED(r_mat_temp,mu_hat,sigma_hat);

p_bar = p_bar_not_norm/sum(p_bar_not_norm);

excess_return = r_mat_temp - mu_hat;
%% Build X and y - OLD
% 
% X = revADm(zeros(n_m+1,length(r_mat_temp)));
% y = revADm(zeros(n_m+1,1));

% X = revADm(zeros(n_m+2,length(r_mat_temp)));
% y = revADm(zeros(n_m+2,1));

% for i_mc=1:n_m 
%     if(1==i_mc)
%         x = p_bar.*excess_return;
%         X(1,:) = x';
%         y(1)   = mc(1)-sum(x);
%     else
%         x = x.*excess_return;
%         
% %         X(i_mc-1,:) = x';
% %         y(i_mc-1)   = mc(i_mc-1) - sum(x);
%         
%         X(i_mc,:) = x';
%         y(i_mc)   = mc(i_mc) - sum(x);
%     end
% end
% tilde_d = (1-exp(r_f)./exp(r_mat_temp));
% d       = tilde_d.*p_bar;   

% y(size(y,1)) = -sum(d);

% X(size(X,1),:) = p_bar';

% X(size(X,1)-1,:) = p_bar';
% X(size(X,1),:)   = d';

%% WITH d AND NOT FIRST MOMENT
% X = forADm(zeros(n_m+1,length(r_mat_temp)),zeros(n_m+1,length(r_mat_temp)));
% y = forADm(zeros(n_m+1,1),zeros(n_m+1,1));

X = forADm(zeros(2,length(r_mat_temp)),zeros(2,length(r_mat_temp)));
y = forADm(zeros(2,1),zeros(2,1));


tilde_d = (1 - exp(r_f*sigma_scale^2)./exp(r_mat_temp) );
d       = tilde_d.*p_bar; 

y(1) = -sum(d);
X(1,:) = d';

% x = p_bar.*excess_return;
% for i_mc=2:n_m 
%         x = x.*excess_return;
%         X(i_mc,:) = x'/((sigma_hat.val)^i_mc);
%         y(i_mc)   = (mc(i_mc) - sum(x))/((sigma_hat.val)^i_mc);
% end
X(size(X,1),:) = p_bar';
%% Solve

epsilon = f_epsilon(X,y);
% epsilon = X'*((X*X')\y);

v1 = (1+epsilon);
v2 = (p_bar.*exp(-r_mat_temp));

q = (v1.*v2)/(v1'*v2);
q_exp = 0;
if(min(q.val)<0)
    q_exp = 1;
    q_exceptions_neg = q_exceptions_neg + 1;
    v2 = (p_bar.*exp(-r_mat_temp));
    v1 = ones(length(p_bar.val),1);
    q = (v1.*v2)/(v1'*v2);
    if(sum(isnan(q.val))>0 || sum(isnan(q.der))>0)
        error('AJ f_node_q')
    end
end

if(sum(isnan(q.val))>0 || sum(isnan(q.der))>0)
    q_exp = 2;
    q_exceptions_nan = q_exceptions_nan + 1;
    v2 = (p_bar.*exp(-r_mat_temp));
    v1 = ones(length(p_bar.val),1);
    q = (v1.*v2)/(v1'*v2);
    if(sum(isnan(q.val))>0 || sum(isnan(q.der))>0)
        error('AJ f_node_q')
    end
end

end
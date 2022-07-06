function [Qd, Qd_exp] = f_Qd(FROM,TO,I_moat,i_ini,r_mat,r_f,mu_hat,sigma_hat,n_m,sigma_id,sigma_scale)


Qd = revADm(zeros(size(FROM,1),length(r_f)));
Qd_exp = zeros(size(FROM,1),length(r_f));
Qd(i_ini,1) = revADm(1);
qd = Qd(:,1);
out_iter = length(r_f);
for t=1:length(r_f)-1
    [Qt,Qt_exp] = f_Qt(FROM(:,t),TO(:,t),t,I_moat,r_mat,r_f(t),mu_hat(:,t),sigma_hat(:,t),n_m,sigma_id,sigma_scale);
    qdt = Qt*qd;
    
    Qd(:,t+1) = qdt;
    Qd_exp(:,t) = Qt_exp;
    qd=qdt;
    out_iter = out_iter -1
%     toc
end


end
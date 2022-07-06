function [Q_t,Q_exp] = f_Qt(FROM_col,TO_col,t,I_moat,r_mat,r_f,mu_hat,sigma_hat,n_m,sigma_id,sigma_scale)

global tape
global tape_cur_id
tic

tape_cur_id_enter = tape_cur_id;

num_independent_var = numel(sigma_id);

Q_t     = revADm(zeros(size(FROM_col,1),size(FROM_col,1)));
Q_exp   = zeros(size(FROM_col,1),1);

cell_T = cell(length(FROM_col),1);
iter = 1;
for i_row=1:length(FROM_col);
    if(sum(i_row == I_moat)==1)
        Q_t(i_row,i_row) = revADm(1,0);
    else
        from = FROM_col(i_row,1);
        to   = TO_col(i_row,1)  ;
        
        [q_t_for,q_exp] = f_node_q(from,to,r_mat(i_row,:),r_f,mu_hat(i_row),sigma_hat(i_row),n_m,sigma_scale(i_row,t));
        
        if(q_exp>0)
            Q_exp(i_row,1) = q_exp;
        end
        
        cur_sigma_id = sigma_id(i_row,t);
        T  = sparse(length(q_t_for.val),num_independent_var);
        T(:,cur_sigma_id) = q_t_for.der;
        cell_T{iter,1} = T;
        iter = iter + 1;
        
        len_q = length(q_t_for.val);
        
        Q_t(from:to,i_row) = revADm(q_t_for.val,tape_cur_id_enter + (1:1:len_q)');
        tape_cur_id_enter = tape_cur_id_enter + len_q;
    end
end

T_send = cell2mat(cell_T);
extendtape(T_send,'ext_Qt')
end



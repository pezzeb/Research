function [gc,gp] = interpol_pricing(Oc,Op,Qd_in,StrikeL,T,rf,AD_scalar)
%INTERPOL_PRICING Price call and put options with interpolation
%   Oc      is the set of call options 
%   Op      is the set of put  options 
%   T       is the grid times
%   Disc    is the discount factors for all grid times
%   StrikL  is the strike levels in the grid
%   Qd      is the complete set of all riskneutral distributions

global tape_cur_id;

enter_length = tape_cur_id;

not_sorted_O = [[Oc,ones(size(Oc,1),1)];[Op,zeros(size(Op,1),1)]];
O = sortrows(not_sorted_O,[2,1]);

delta_T = f_Delta_T(T);
DiscT = [exp(-delta_T.*rf(1:end-1));1];

O(:,2) = round(1000000*O(:,2))/1000000;

u_t_time_finder = unique(O(:,2));

[~,used_time_id] = ismember(u_t_time_finder,round(1000000*T)/1000000);
Disc = DiscT(used_time_id);

Qd = Qd_in(:,used_time_id);

%iteration setup
call_id = 1;
put_id  = 1;
iter    = 1;
nO      = size(O,1);
nOc     = sum(O(:,4));
nOp     = nO-nOc;



if(strcmp(AD_scalar,'AD'))
    gc = revADm(zeros(nOc,1));
    gp = revADm(zeros(nOp,1));
elseif(strcmp(AD_scalar,'scalar'))
    gc = zeros(nOc,1);
    gp = zeros(nOp,1);
end


cell_T = {};

for t=1:length(u_t_time_finder)
    enter_time = u_t_time_finder(t);
    while(1+iter <= nO && O(iter,2) == enter_time)
        fprintf('iteration in Interpolation pricing: %0.f \n',size(O,1)-iter)
        qd      = Qd(:,t);
        Disc_t  = Disc(t);
        Ot      = O(iter,1);
        
        upper_idx  = find(StrikeL>Ot,1,'first');
        lower_idx  = upper_idx - 1;
        
        Kl = StrikeL(lower_idx);
        Ku = StrikeL(upper_idx);
        
        c_u        = max(StrikeL-Ku,0);
        c_l        = max(StrikeL-Kl,0);
        p_u        = max(Ku-StrikeL,0);
        p_l        = max(Kl-StrikeL,0);
        
        Cupper_val = c_u'*qd.val*Disc_t;
        Clower_val = c_l'*qd.val*Disc_t;
        
        Pupper_val = p_u'*qd.val*Disc_t;
        Plower_val = p_l'*qd.val*Disc_t;
        
        if(O(iter,1)==O(iter+1,1))
            %Put & Call

            [gct,gpt,T] = help_func(Cupper_val,Clower_val,Pupper_val,Plower_val,Disc_t,StrikeL,qd,Ot,'cp',c_u,c_l,p_u,p_l,enter_length,AD_scalar);
            gc(call_id,1) = gct;
            gp(put_id,1)  = gpt;
            
            call_id = call_id + 1;
            put_id  = put_id  + 1;
            
            cell_T{end+1,1} = T;
            
            iter = iter + 2;
        elseif(1==O(iter,4))
            %CALL
            
            [gct,~,T] = help_func(Cupper_val,Clower_val,-1,-1,Disc(t),StrikeL,Qd(:,t),Ot,'c',c_u,c_l,p_u,p_l,enter_length,AD_scalar);
            gc(call_id,1) = gct;
            
            call_id = call_id + 1;
            
            cell_T{end+1,1} = T;
            
            iter = iter + 1;
        elseif(0==O(iter,4))
            %PUT
            
            [~,gpt,T] = help_func(-1,-1,Pupper_val,Plower_val,Disc_t,StrikeL,qd,Ot,'p',c_u,c_l,p_u,p_l,enter_length,AD_scalar);
            gp(put_id,1) = gpt;
            
            put_id  = put_id  + 1;
            
            cell_T{end+1,1} = T;
            
            iter = iter + 1;
        else 
            error('ERROR in Interpol pricing');
        end
    end
    
    if(iter<=nO)
        qd      = Qd(:,t);
        Disc_t  = Disc(t);
        Ot      = O(iter,1);
        
        upper_idx  = find(StrikeL>Ot,1,'first');
        lower_idx  = upper_idx - 1;
        
        Kl = StrikeL(lower_idx);
        Ku = StrikeL(upper_idx);
        
        c_u        = max(StrikeL-Ku,0);
        c_l        = max(StrikeL-Kl,0);
        p_u        = max(Ku-StrikeL,0);
        p_l        = max(Kl-StrikeL,0);
        
        Cupper_val = c_u'*qd.val*Disc_t;
        Clower_val = c_l'*qd.val*Disc_t;
        
        Pupper_val = p_u'*qd.val*Disc_t;
        Plower_val = p_l'*qd.val*Disc_t;
        
        if(iter == nO)
            if(1==O(iter,4))
                %CALL
                
                [gct,~] = help_func(Cupper_val,Clower_val,-1,-1,Disc(t),StrikeL,Qd(:,t),Ot,'c',c_u,c_l,p_u,p_l,enter_length,AD_scalar);
                gc(call_id,1) = gct;
                
                call_id = call_id + 1;
                
                cell_T{end+1,1} = T;
                
                iter = iter + 1;
            elseif(0==O(iter,4))
                %PUT
                
                [~,gpt] = help_func(-1,-1,Pupper_val,Plower_val,Disc_t,StrikeL,qd,Ot,'p',c_u,c_l,p_u,p_l,enter_length,AD_scalar);
                gp(put_id,1) = gpt;
                
                put_id  = put_id  + 1;
                
                cell_T{end+1,1} = T;
                
                iter = iter + 1;
            else
                error('ERROR in Interpol pricing');
            end
        end
    else
    end
end

T_send = cell2mat(cell_T);
if(strcmp(AD_scalar,'AD'))
    extendtape(T_send,'interpol_pricing');
elseif(strcmp(AD_scalar,'scalar'))
    %nothing
end
end

function [gc,gp,T] = help_func(Cupper_val,Clower_val,Pupper_val,Plower_val,Disc_t,grid,qd,K,cp_case,c_u,c_l,p_u,p_l,enter_length,AD_scalar)
global tape_cur_id

upper_idx  = find(grid>K,1,'first');
lower_idx  = upper_idx - 1;

Kl = grid(lower_idx);
Ku = grid(upper_idx);

d_lower = K-Kl;
d_upper = K-Ku;

qd_val = qd.val;
short_sum = sum(qd_val(upper_idx:length(qd_val)));  %CANNOT USE end, not implemented in RevADm
long_sum  = short_sum + qd_val(lower_idx);

N_lower = 2/(grid(lower_idx+1)-grid(lower_idx-1));
N_upper = 2/(grid(upper_idx+1)-grid(upper_idx-1));

dC2_lower_val = N_lower*qd.val(lower_idx);
dC2_upper_val = N_upper*qd.val(upper_idx);

lambda = (Ku-K)/(Ku-Kl);

if(strcmp(cp_case,'cp'))
    
    T = sparse(2,enter_length);
    
    dC_lower_val = -Disc_t*long_sum;
    dC_upper_val = -Disc_t*short_sum;
    
    dP_lower = Disc_t*(1-long_sum);
    dP_upper = Disc_t*(1-short_sum);
    
    gc_val = lambda*(Clower_val + dC_lower_val*d_lower + 1/2*dC2_lower_val*d_lower.^2) + ...
        (1-lambda)*(Cupper_val  + dC_upper_val*d_upper + 1/2*dC2_upper_val*d_upper.^2);
    
    gp_val =  lambda*(Plower_val + dP_lower*d_lower + 1/2*dC2_lower_val*d_lower.^2) + ...
        (1-lambda)*(Pupper_val   + dP_upper*d_upper + 1/2*dC2_upper_val*d_upper.^2);
    
    q_id = qd.id;
    
    T(1,q_id') = Disc_t*(c_l*lambda + (1-lambda)*c_u);
    T(2,q_id') = Disc_t*(p_l*lambda + (1-lambda)*p_u);
    
    T([1;2],[q_id(upper_idx:end)',q_id(lower_idx)])   = T([1;2],[q_id(upper_idx:end)',q_id(lower_idx)]) - lambda*Disc_t*d_lower;
    T([1;2],[q_id(upper_idx:end)'])                   = T([1;2],[q_id(upper_idx:end)'])                 - (1-lambda)*Disc_t*d_upper;
    T([1;2],[q_id(lower_idx),q_id(upper_idx)])        = T([1;2],[q_id(lower_idx),q_id(upper_idx)])      + [N_lower*lambda/2*d_lower.^2,N_upper*(1-lambda)/2*d_upper.^2;N_lower*lambda/2*d_lower.^2,N_upper*(1-lambda)/2*d_upper.^2];
    
    if(strcmp(AD_scalar,'AD'))
        gc = revADm(gc_val,tape_cur_id+1);
        gp = revADm(gp_val,tape_cur_id+2);
        tape_cur_id = tape_cur_id + 2;
    elseif(strcmp(AD_scalar,'scalar'))
        gc = gc_val;
        gp = gp_val;
    end
    
    
elseif(strcmp(cp_case,'c'))
    T = sparse(1,enter_length);
    dC_lower_val = -Disc_t*long_sum;
    dC_upper_val = -Disc_t*short_sum;
    
    gp = -1;
    
    gc_val =   lambda*(Clower_val + dC_lower_val*d_lower + 1/2*dC2_lower_val*d_lower.^2) + ...
        (1-lambda)*(Cupper_val    + dC_upper_val*d_upper + 1/2*dC2_upper_val*d_upper.^2);
    
    if(strcmp(AD_scalar,'AD'))
        gc = revADm(gc_val,tape_cur_id+1);
        tape_cur_id = tape_cur_id + 1;
    elseif(strcmp(AD_scalar,'scalar'))
        gc = gc_val;
    end
    
    q_id = qd.id;
    
    T(1,q_id') = Disc_t*(c_l*lambda + (1-lambda)*c_u);
    
    T(1,[q_id(upper_idx:end)',q_id(lower_idx)])   = T(1,[q_id(upper_idx:end)',q_id(lower_idx)]) - [lambda*Disc_t*d_lower];
    T(1,[q_id(upper_idx:end)'])                   = T(1,[q_id(upper_idx:end)'])                 - [(1-lambda)*Disc_t*d_upper];
    T(1,[q_id(lower_idx),q_id(upper_idx)])        = T(1,[q_id(lower_idx),q_id(upper_idx)])      + [N_lower*lambda/2*d_lower.^2,N_upper*(1-lambda)/2*d_upper.^2];
    
elseif(strcmp(cp_case,'p'))
    T = sparse(1,enter_length);
    
    dP_lower = Disc_t*(1-long_sum);
    dP_upper = Disc_t*(1-short_sum);
    
    gp_val =  lambda*(Plower_val + dP_lower*d_lower + 1/2*dC2_lower_val*d_lower.^2) + ...
        (1-lambda)*(Pupper_val   + dP_upper*d_upper + 1/2*dC2_upper_val*d_upper.^2);
    
    if(strcmp(AD_scalar,'AD'))
        gp = revADm(gp_val,tape_cur_id+1);
        tape_cur_id = tape_cur_id + 1;
    elseif(strcmp(AD_scalar,'scalar'))
        gp = gp_val;
    end
    
    gc = -1;
    
    q_id = qd.id;
    
    T(1,q_id') = Disc_t*(p_l*lambda + (1-lambda)*p_u);
    
    T(1,[q_id(upper_idx:end)',q_id(lower_idx)])   = T(1,[q_id(upper_idx:end)',q_id(lower_idx)]) - [lambda*Disc_t*d_lower];
    T(1,[q_id(upper_idx:end)'])                   = T(1,[q_id(upper_idx:end)'])                 - [(1-lambda)*Disc_t*d_upper];
    T(1,[q_id(lower_idx),q_id(upper_idx)])        = T(1,[q_id(lower_idx),q_id(upper_idx)])      + [N_lower*lambda/2*d_lower.^2,N_upper*(1-lambda)/2*d_upper.^2];
    
end
end
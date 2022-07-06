function [ E_out ] = f_EF(EF_in, Oc,Op,type,start_vol,Sini,EF_con,Tc)
%F_E Returns a matrix E.
%   PARAMETERS:
%   EF      - is a user defined E matrix that is choosen as default
%   Oc      - is the call option object
%   Op      - is the put  option object
%   type    - is the type of E

Ok = [Oc(:,1);Op(:,1)];
Ot = [Oc(:,2);Op(:,2)];
K  = unique(Ok);
T  = unique(Ot);

if(strcmp(type,'constant'))
    n = size(Oc,1) + size(Op,1);
    v_EF_in = ones(n,1).*EF_in;
    E_out = spdiags(v_EF_in,0,sparse(n,n));
elseif(strcmp(type,'normal_dist'))
    %     unique_maturity = unique(Oc(:,2));
    E_diag =[];
    %     vol_scale = start_vol*sqrt(unique_maturity);
    r_t = normpdf(log(K/Sini),0,start_vol);
    r_scale = r_t/min(r_t);
    
    [~,idx] = ismember(Ok,K);
    r_out = r_scale(idx);
    
    E = EF_con*exp(25*(r_out-1));
    E_out = sparse(diag(E));
    
elseif(strcmp(type,'K'))
%     for k=1:length(K)
%         n_K(k,1) = sum(K(k)==Ok);
%     end
%     n_tot = sum(n_K);
%     ratio_K = n_K/n_tot;
%     
%         for j=1:length(Ok)
%             ratio_K_scale(j,1) = 1/ratio_K(K==Ok(j)); %1/(n_K(K==Ok(j))^2);
%         end
%         
%     ratio_K_scale = ratio_K_scale / max(ratio_K_scale);
%     E_scale = ratio_K_scale.*EF_con;
%     
%     E_out = sparse(diag(E_scale));
%     

K_max = max(K);
K_min = min(K);

K_disc = K_min:(K_max-K_min)/20:K_max;
E = zeros(size(Ok));
for i=1:length(K_disc)
    if(length(K_disc)==i)
        n       = numel(Ok)-sum(n);
        ratio   = n/numel(Ok);
        E(logical(K_disc(i) == Ok)) = ratio;
    else
        number_id   = (K_disc(i) <= Ok).*(Ok < K_disc(i+1));
        n           = sum(number_id);
        ratio       = n/numel(Ok);
        E(logical((K_disc(i) <= Ok).*(Ok < K_disc(i+1)))) = ratio;
    end
    
end
E_scale = E/max(E).*EF_con;
E_out = sparse(diag(E_scale));

elseif(strcmp(type,'T'))
    
    [~,idx] = ismember(T,Tc);
    
    for k=1:length(Ot)
        E_temp(k,1) = idx((Ot(k)==T))/length(Tc);
    end
    E_scale = E_temp/min(E_temp)*EF_con;
    E_out = sparse(diag(E_scale));
elseif(strcmp(type,'K and T'))
    
    for k=1:length(K)
        n_K(k,1) = sum(K(k)==Ok);
    end
    n_tot = sum(n_K);
    ratio_K = n_K/n_tot;
    
    for j=1:length(Ok)
        E_tempT(j,1) = 1/(n_K(K==Ok(j))^2);
    end
    
    [~,idx] = ismember(T,Tc);
    
    for k=1:length(Ot)
        E_tempK(k,1) = idx((Ot(k)==T))/length(Tc);
    end
    
    E_temp = E_tempT.*E_tempK;
    E_scale = E_temp/min(E_temp)*EF_con;
    E_out = sparse(diag(E_scale));
end
end




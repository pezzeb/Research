function [Hh_output,h11,h12,h13,h14,h15] = f_Hh(T,K,gamma_m,phi_m,xi_m,upsilon_m,theta_m,Hh_type,u_bar_val)
%F_Hh Returns the matrix Hh. 
%   PARAMETER:
%   delta_T - (vector) is a increments of the grids maturity times
%   delta_K - (vector) is a increments of the grids strike prices
%
%   gamma   - (matrix) is the gamma matrix of the problem
%   phi     - (matrix) is the phi matrix of the problem
%   xi      - (matrix) is the xi matrix of the problem
%   upsilon - (matrix) is the upsilon matrix of the problem
%   theta   - (matrix) is the theta matrix of the problem
%
%   RETURN: 
%   Hh      - (matrix) is the Hh matrix.

delta_T = f_Delta_T(T);
delta_K = diff(K);

nT = length(delta_T);
nK = length(delta_K);

gamma   = vec_v(gamma_m);
phi     = vec_v(phi_m);
xi      = xi_m;
upsilon = vec_v(upsilon_m);


%% create hat
gamma_hat  = gamma./repmat(delta_T,size(gamma_m,1),1);
phi_hat    = phi./repmat(delta_K,size(phi_m,2),1);
xi_hat     = vec_h(2*xi./repmat((delta_T(2:end).^2 .* delta_T(1:end-1).^2 .* (delta_T(2:end) + delta_T(1:end-1)) )',size(xi_m,1),1));
upsilon_hat= 2*upsilon./repmat((delta_K(2:end).^2.*delta_K(1:end-1).^2.*(delta_K(2:end) + delta_K(1:end-1))),size(upsilon_m,2),1);

theta_m  = theta_m./((delta_K(2:end) + delta_K(1:end-1))*(delta_T(2:end) + delta_T(1:end-1))');
theta_hat   = vec_v(theta_m);
%% create B
P = f_P(nK,nT);

B_gamma = f_B(-ones(1,nT),ones(1,nT));
B_phi   = f_B(-ones(1,nK),ones(1,nK));
B_xi    = f_B(delta_T(2:end),-(delta_T(2:end)+delta_T(1:end-1)),delta_T(1:end-1));
B_upsilon = f_B(delta_K(2:end),-(delta_K(2:end)+delta_K(1:end-1)),delta_K(1:end-1));
B_theta = f_B(-ones(1,nK-1),zeros(1,nK-1),ones(1,nK-1));
%% create A

A_gamma   = sparse(f_A(B_gamma,nK+1));
A_phi     = sparse(f_A(B_phi,nT+1));
A_xi      = sparse(f_A(B_xi,nK+1));
A_upsilon = sparse(f_A(B_upsilon,nT+1));

A_theta_temp = f_A(B_theta,nT-1);
A_theta = [sparse(size(A_theta_temp,1),2*(nK+1)),A_theta_temp] - [A_theta_temp,sparse(size(A_theta_temp,1),2*(nK+1))];

if(strcmp(Hh_type,'different'))
%     h1 = P'*A_gamma'*diag(gamma_hat(:))*A_gamma*P;
%     h2 = A_phi'*diag(phi_hat(:))*A_phi;
%     h3 = P'*A_xi'*diag(xi_hat(:))*A_xi*P;
%     h4 = A_upsilon'*diag(upsilon_hat(:))*A_upsilon;
%     h5 = A_theta'*diag(theta_hat(:))*A_theta;
    
    len11= length(gamma_hat);
    D11  = spdiags(gamma_hat(:),0,sparse(len11,len11));
    h11  = P'*A_gamma'*D11*A_gamma*P;
    
    len12= length(phi_hat);
    D12 = spdiags(phi_hat(:),0,sparse(len12,len12));
    h12 = A_phi'*D12*A_phi;
    
    len13= length(xi_hat);
    D13 = spdiags(xi_hat(:),0,sparse(len13,len13));
% D13 = 1;
     
    h13 = P'*(A_xi'*D13*A_xi)*P;
    
    len14= length(upsilon_hat);
    D14 = spdiags(upsilon_hat(:),0,sparse(len14,len14));
    h14 = A_upsilon'*D14*A_upsilon;
    
    len15= length(theta_hat);
    D15 = spdiags(theta_hat(:),0,sparse(len15,len15));
    h15 = A_theta'*D15*A_theta;
    
elseif(strcmp(Hh_type,'constant'))
    APg = A_gamma*P;
    APx = A_xi*P;
    
    h1 = APg'*APg*gamma_hat(1,1);
    fprintf('The matrix h1  WAS COMPLETED in: %f seconds\n',toc);
    h2 = A_phi'*A_phi*phi_hat(1,1);
    fprintf('The matrix h2  WAS COMPLETED in: %f seconds\n',toc);
    h3 = APx'*APx*xi_hat(1,1);
    fprintf('The matrix h3  WAS COMPLETED in: %f seconds\n',toc);
    h4 = A_upsilon'*A_upsilon*upsilon_hat(1,1);
    fprintf('The matrix h4  WAS COMPLETED in: %f seconds\n',toc);
    h5 = A_theta'*A_theta*theta_hat(1,1);
    fprintf('The matrix h5  WAS COMPLETED in: %f seconds\n',toc);
    
else
    error('No such type i f_Hh');
end
%     Hh_output = h1+h2+h3+h4+h5;
    Hh_output = h11+h12+h13+h14+h15;


end

function [P_out] = f_P(nK,nT)

P_out = sparse(nK+1,nT+1);

for i=1:(nK+1)*(nT+1) 
    j = 1 + floor((i-1)/(nT+1)) + (i-1)*(nK+1)-(nK+1)*(nT+1)*floor((i-1)/(nT+1));
    P_out(i,j) = 1;
end
end

function [B_out] = f_B(diagonal,super_diagonal,super_super_diagonal)

%ERROR
%first or second term
if nargin == 2 && (length(diagonal) == length(super_diagonal))
    B_out = zeros(length(diagonal),length(diagonal)+1);
    B_out(:,1:end-1) = B_out(:,1:end-1) + diag(diagonal);
    B_out(:,2:end)   = B_out(:,2:end  ) + diag(super_diagonal);
elseif nargin == 3 && (length(diagonal) == length(super_diagonal) && length(diagonal) == length(super_super_diagonal))
    B_out = zeros(length(diagonal),length(diagonal)+2);
    B_out(:,1:end-2) = B_out(:,1:end-2) + diag(diagonal);
    B_out(:,2:end-1) = B_out(:,2:end-1) + diag(super_diagonal);
    B_out(:,3:end  ) = B_out(:,3:end  ) + diag(super_super_diagonal);
else
    %ERROR
end
end

function [A_out] = f_A(B,n_rep)

A_out = zeros(size(B)*n_rep);

h = size(B,1);
w = size(B,2);

for i=0:n_rep-1
    A_out(1 + i*h: h + i*h, 1+ i*w: w+i*w) = B;
end
end
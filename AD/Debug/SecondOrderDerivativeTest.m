%Add the pricing functionality folder

x = forADm(5,[1,0]);
y = forADm(4,[0,1]);

% SSE = x^2 * y^2;

%%
x_hes = forADm(x,[1 0]);
xsvar = x_hes^2*y^2;

xsvar.val
xder_svar = xsvar.der
xder_svar.der

%%
y_hes = forADm(y,[0 1]);
ysvar = y_hes^2*x^2;

ysvar.val
yder_svar = ysvar.der
yder_svar.der

%%
S0_  = 100;
K_   = 100;
T_   = 1;
r_   = 0;
div_ = 0;

nu0_    = 0.01;
kappa_  = 3;
eta_    = 0.01;
theta_  = 0.3;
rho_    = -0.2;
etabar_ = kappa_*eta_;
Vtheta_ = theta_^2;

nu0     = forADm(nu0_     ,[1 0 0 0 0 ]);
kappa   = forADm(kappa_   ,[0 1 0 0 0 ]);
etabar  = forADm(etabar_  ,[0 0 1 0 0 ]);
Vtheta  = forADm(Vtheta_  ,[0 0 0 1 0 ]);
rho     = forADm(rho_     ,[0 0 0 0 1 ]);

nu0_hes    = forADm(nu0   ,[1 0 0 0 0 ]);
kappa_hes  = forADm(kappa ,[0 1 0 0 0 ]);
etabar_hes = forADm(etabar,[0 0 1 0 0 ]);
Vtheta_hes = forADm(Vtheta,[0 0 0 1 0 ]);
rho_hes    = forADm(rho   ,[0 0 0 0 1 ]);

grad = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0    ;kappa;etabar;Vtheta;rho},S0_,K_,T_,r_,div_,[])
%% lite lättare?
hes_e1 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_hes;kappa_;etabar_;Vtheta_;rho_},S0_,K_,T_,r_,div_,[])
hes_e2 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_hes;etabar_;Vtheta_;rho_},S0_,K_,T_,r_,div_,[])
hes_e3 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_hes;Vtheta_;rho_},S0_,K_,T_,r_,div_,[])
hes_e4 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_;Vtheta_hes;rho_},S0_,K_,T_,r_,div_,[])
hes_e5 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_;Vtheta_;rho_hes},S0_,K_,T_,r_,div_,[])
%% svår
hes_h1 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_hes;kappa;etabar;Vtheta;rho},S0_,K_,T_,r_,div_,[])
hes_h2 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0;kappa_hes;etabar;Vtheta;rho},S0_,K_,T_,r_,div_,[])
hes_h3 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0;kappa;etabar_hes;Vtheta;rho},S0_,K_,T_,r_,div_,[])
hes_h4 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0;kappa;etabar;Vtheta_hes;rho},S0_,K_,T_,r_,div_,[])
hes_h5 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0;kappa;etabar;Vtheta;rho_hes},S0_,K_,T_,r_,div_,[])
%% jättesvår
hes_vh   = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_hes;kappa_hes;etabar_hes;Vtheta_hes;rho_hes},S0_,K_,T_,r_,div_,[])

%% Numierska derivator
%nu0 second order
dnu0 = 0.00001;
dkap = 0.00001;
deba = 0.00001;
dthe = 0.00001;
drho = 0.00001;

fm_    = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);

fu_nu0 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_+dnu0;kappa_;etabar_;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fd_nu0 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_-dnu0;kappa_;etabar_;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fu_kap = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_+dkap;etabar_;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fd_kap = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_-dkap;etabar_;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fu_eba = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_+deba;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fd_eba = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_-deba;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fu_the = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_;Vtheta_+dthe;rho_},S0_,K_,T_,r_,div_,[]);
fd_the = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_;Vtheta_-dthe;rho_},S0_,K_,T_,r_,div_,[]);
fu_rho = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_;Vtheta_;rho_+drho},S0_,K_,T_,r_,div_,[]);
fd_rho = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_;Vtheta_;rho_-drho},S0_,K_,T_,r_,div_,[]);
%% Nu0 - Second order
fu_nu0_u_kap = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_+dnu0;kappa_+dkap;etabar_;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fu_nu0_u_eba = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_+dnu0;kappa_;etabar_+deba;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fu_nu0_u_the = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_+dnu0;kappa_;etabar_;Vtheta_+dthe;rho_},S0_,K_,T_,r_,div_,[]);
fu_nu0_u_rho = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_+dnu0;kappa_;etabar_;Vtheta_;rho_+drho},S0_,K_,T_,r_,div_,[]);

fu_nu0_d_kap = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_+dnu0;kappa_-dkap;etabar_;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fu_nu0_d_eba = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_+dnu0;kappa_;etabar_-deba;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fu_nu0_d_the = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_+dnu0;kappa_;etabar_;Vtheta_-dthe;rho_},S0_,K_,T_,r_,div_,[]);
fu_nu0_d_rho = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_+dnu0;kappa_;etabar_;Vtheta_;rho_-drho},S0_,K_,T_,r_,div_,[]);

fd_nu0_u_kap = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_-dnu0;kappa_+dkap;etabar_;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fd_nu0_u_eba = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_-dnu0;kappa_;etabar_+deba;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fd_nu0_u_the = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_-dnu0;kappa_;etabar_;Vtheta_+dthe;rho_},S0_,K_,T_,r_,div_,[]);
fd_nu0_u_rho = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_-dnu0;kappa_;etabar_;Vtheta_;rho_+drho},S0_,K_,T_,r_,div_,[]);

fd_nu0_d_kap = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_-dnu0;kappa_-dkap;etabar_;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fd_nu0_d_eba = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_-dnu0;kappa_;etabar_-deba;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fd_nu0_d_the = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_-dnu0;kappa_;etabar_;Vtheta_-dthe;rho_},S0_,K_,T_,r_,div_,[]);
fd_nu0_d_rho = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_-dnu0;kappa_;etabar_;Vtheta_;rho_-drho},S0_,K_,T_,r_,div_,[]);
%% Kappa - Second order
fu_kap_u_nu0 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_+dnu0;kappa_+dkap;etabar_;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fu_kap_u_eba = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_+dkap;etabar_+deba;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fu_kap_u_the = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_+dkap;etabar_;Vtheta_+dthe;rho_},S0_,K_,T_,r_,div_,[]);
fu_kap_u_rho = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_+dkap;etabar_;Vtheta_;rho_+drho},S0_,K_,T_,r_,div_,[]);

fu_kap_d_nu0 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_-dnu0;kappa_+dkap;etabar_;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fu_kap_d_eba = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_+dkap;etabar_-deba;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fu_kap_d_the = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_+dkap;etabar_;Vtheta_-dthe;rho_},S0_,K_,T_,r_,div_,[]);
fu_kap_d_rho = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_+dkap;etabar_;Vtheta_;rho_-drho},S0_,K_,T_,r_,div_,[]);

fd_kap_u_nu0 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_+dnu0;kappa_-dkap;etabar_;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fd_kap_u_eba = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_-dkap;etabar_+deba;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fd_kap_u_the = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_-dkap;etabar_;Vtheta_+dthe;rho_},S0_,K_,T_,r_,div_,[]);
fd_kap_u_rho = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_-dkap;etabar_;Vtheta_;rho_+drho},S0_,K_,T_,r_,div_,[]);

fd_kap_d_nu0 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_-dnu0;kappa_-dkap;etabar_;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fd_kap_d_eba = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_-dkap;etabar_-deba;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fd_kap_d_the = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_-dkap;etabar_;Vtheta_-dthe;rho_},S0_,K_,T_,r_,div_,[]);
fd_kap_d_rho = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_-dkap;etabar_;Vtheta_;rho_-drho},S0_,K_,T_,r_,div_,[]);
%% etabar - Second order
fu_eba_u_nu0 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_+dnu0;kappa_;etabar_+deba;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fu_eba_u_kap = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_+dkap;etabar_+deba;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fu_eba_u_the = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_+deba;Vtheta_+dthe;rho_},S0_,K_,T_,r_,div_,[]);
fu_eba_u_rho = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_+deba;Vtheta_;rho_+drho},S0_,K_,T_,r_,div_,[]);

fu_eba_d_nu0 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_-dnu0;kappa_;etabar_+deba;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fu_eba_d_kap = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_-dkap;etabar_+deba;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fu_eba_d_the = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_+deba;Vtheta_-dthe;rho_},S0_,K_,T_,r_,div_,[]);
fu_eba_d_rho = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_+deba;Vtheta_;rho_-drho},S0_,K_,T_,r_,div_,[]);

fd_eba_u_nu0 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_+dnu0;kappa_;etabar_-deba;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fd_eba_u_kap = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_+dkap;etabar_-deba;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fd_eba_u_the = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_-deba;Vtheta_+dthe;rho_},S0_,K_,T_,r_,div_,[]);
fd_eba_u_rho = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_-deba;Vtheta_;rho_+drho},S0_,K_,T_,r_,div_,[]);

fd_eba_d_nu0 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_-dnu0;kappa_;etabar_-deba;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fd_eba_d_kap = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_-dkap;etabar_-deba;Vtheta_;rho_},S0_,K_,T_,r_,div_,[]);
fd_eba_d_the = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_-deba;Vtheta_-dthe;rho_},S0_,K_,T_,r_,div_,[]);
fd_eba_d_rho = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_-deba;Vtheta_;rho_-drho},S0_,K_,T_,r_,div_,[]);
%% Vtheta
fu_the_u_nu0 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_+dnu0;kappa_;etabar_;Vtheta_+dthe;rho_},S0_,K_,T_,r_,div_,[]);
fu_the_u_kap = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_+dkap;etabar_;Vtheta_+dthe;rho_},S0_,K_,T_,r_,div_,[]);
fu_the_u_eba = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_+deba;Vtheta_+dthe;rho_},S0_,K_,T_,r_,div_,[]);
fu_the_u_rho = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_;Vtheta_+dthe;rho_+drho},S0_,K_,T_,r_,div_,[]);

fu_the_d_nu0 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_-dnu0;kappa_;etabar_;Vtheta_+dthe;rho_},S0_,K_,T_,r_,div_,[]);
fu_the_d_kap = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_-dkap;etabar_;Vtheta_+dthe;rho_},S0_,K_,T_,r_,div_,[]);
fu_the_d_eba = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_-deba;Vtheta_+dthe;rho_},S0_,K_,T_,r_,div_,[]);
fu_the_d_rho = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_;Vtheta_+dthe;rho_-drho},S0_,K_,T_,r_,div_,[]);

fd_the_u_nu0 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_+dnu0;kappa_;etabar_;Vtheta_-dthe;rho_},S0_,K_,T_,r_,div_,[]);
fd_the_u_kap = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_+dkap;etabar_;Vtheta_-dthe;rho_},S0_,K_,T_,r_,div_,[]);
fd_the_u_eba = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_+deba;Vtheta_-dthe;rho_},S0_,K_,T_,r_,div_,[]);
fd_the_u_rho = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_;Vtheta_-dthe;rho_+drho},S0_,K_,T_,r_,div_,[]);

fd_the_d_nu0 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_-dnu0;kappa_;etabar_;Vtheta_-dthe;rho_},S0_,K_,T_,r_,div_,[]);
fd_the_d_kap = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_-dkap;etabar_;Vtheta_-dthe;rho_},S0_,K_,T_,r_,div_,[]);
fd_the_d_eba = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_-deba;Vtheta_-dthe;rho_},S0_,K_,T_,r_,div_,[]);
fd_the_d_rho = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_;Vtheta_-dthe;rho_-drho},S0_,K_,T_,r_,div_,[]);
%% Rho
fu_rho_u_nu0 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_+drho;kappa_;etabar_;Vtheta_;rho_+drho},S0_,K_,T_,r_,div_,[]);
fu_rho_u_kap = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_+dkap;etabar_;Vtheta_;rho_+drho},S0_,K_,T_,r_,div_,[]);
fu_rho_u_eba = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_+deba;Vtheta_;rho_+drho},S0_,K_,T_,r_,div_,[]);
fu_rho_u_the = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_;Vtheta_+dthe;rho_+drho},S0_,K_,T_,r_,div_,[]);

fu_rho_d_nu0 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_-drho;kappa_;etabar_;Vtheta_;rho_+drho},S0_,K_,T_,r_,div_,[]);
fu_rho_d_kap = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_-dkap;etabar_;Vtheta_;rho_+drho},S0_,K_,T_,r_,div_,[]);
fu_rho_d_eba = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_-deba;Vtheta_;rho_+drho},S0_,K_,T_,r_,div_,[]);
fu_rho_d_the = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_;Vtheta_-dthe;rho_+drho},S0_,K_,T_,r_,div_,[]);

fd_rho_u_nu0 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_+drho;kappa_;etabar_;Vtheta_;rho_-drho},S0_,K_,T_,r_,div_,[]);
fd_rho_u_kap = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_+dkap;etabar_;Vtheta_;rho_-drho},S0_,K_,T_,r_,div_,[]);
fd_rho_u_eba = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_+deba;Vtheta_;rho_-drho},S0_,K_,T_,r_,div_,[]);
fd_rho_u_the = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_;Vtheta_+dthe;rho_-drho},S0_,K_,T_,r_,div_,[]);

fd_rho_d_nu0 = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_-drho;kappa_;etabar_;Vtheta_;rho_-drho},S0_,K_,T_,r_,div_,[]);
fd_rho_d_kap = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_-dkap;etabar_;Vtheta_;rho_-drho},S0_,K_,T_,r_,div_,[]);
fd_rho_d_eba = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_-deba;Vtheta_;rho_-drho},S0_,K_,T_,r_,div_,[]);
fd_rho_d_the = Integralpricing('HestonNewParam','not_alpha_psquad',{nu0_;kappa_;etabar_;Vtheta_-dthe;rho_-drho},S0_,K_,T_,r_,div_,[]);
%% NUMERISKA DERIVATOR
%Nu0
numder1Order(1,1) = (fu_nu0 - fd_nu0)/(2*dnu0);
numder2Order(1,1) = (fu_nu0 - 2*fm_ + fd_nu0)/(dnu0^2);
numder2Order(1,2) = (fu_nu0_u_kap - fu_nu0_d_kap - fd_nu0_u_kap + fd_nu0_d_kap)/(4*dnu0*dkap);
numder2Order(1,3) = (fu_nu0_u_eba - fu_nu0_d_eba - fd_nu0_u_eba + fd_nu0_d_eba)/(4*dnu0*deba);
numder2Order(1,4) = (fu_nu0_u_the - fu_nu0_d_the - fd_nu0_u_the + fd_nu0_d_the)/(4*dnu0*dthe);
numder2Order(1,5) = (fu_nu0_u_rho - fu_nu0_d_rho - fd_nu0_u_rho + fd_nu0_d_rho)/(4*dnu0*drho);
%Kappa
numder1Order(1,2) = (fu_kap - fd_kap)/(2*dkap);
numder2Order(2,2) = (fu_kap - 2*fm_ + fd_kap)/(dkap^2);
numder2Order(2,1) = (fu_kap_u_nu0 - fu_kap_d_nu0 - fd_kap_u_nu0 + fd_kap_d_nu0)/(4*dkap*dnu0);
numder2Order(2,3) = (fu_kap_u_eba - fu_kap_d_eba - fd_kap_u_eba + fd_kap_d_eba)/(4*dkap*deba);
numder2Order(2,4) = (fu_kap_u_the - fu_kap_d_the - fd_kap_u_the + fd_kap_d_the)/(4*dkap*dthe);
numder2Order(2,5) = (fu_kap_u_rho - fu_kap_d_rho - fd_kap_u_rho + fd_kap_d_rho)/(4*dkap*drho);
%etabar
numder1Order(1,3) = (fu_eba - fd_eba)/(2*deba);
numder2Order(3,3) = (fu_eba - 2*fm_ + fd_eba)/(deba^2);
numder2Order(3,1) = (fu_eba_u_nu0 - fu_eba_d_nu0 - fd_eba_u_nu0 + fd_eba_d_nu0)/(4*deba*dnu0);
numder2Order(3,2) = (fu_eba_u_kap - fu_eba_d_kap - fd_eba_u_kap + fd_eba_d_kap)/(4*deba*dkap);
numder2Order(3,4) = (fu_eba_u_the - fu_eba_d_the - fd_eba_u_the + fd_eba_d_the)/(4*deba*dthe);
numder2Order(3,5) = (fu_eba_u_rho - fu_eba_d_rho - fd_eba_u_rho + fd_eba_d_rho)/(4*deba*drho);
%Vtheta
numder1Order(1,4) = (fu_the - fd_the)/(2*dthe);
numder2Order(4,4) = (fu_the - 2*fm_ + fd_the)/(dthe^2);
numder2Order(4,1) = (fu_the_u_nu0 - fu_the_d_nu0 - fd_the_u_nu0 + fd_the_d_nu0)/(4*dthe*dnu0);
numder2Order(4,2) = (fu_the_u_kap - fu_the_d_kap - fd_the_u_kap + fd_the_d_kap)/(4*dthe*dkap);
numder2Order(4,3) = (fu_the_u_eba - fu_the_d_eba - fd_the_u_eba + fd_the_d_eba)/(4*dthe*deba);
numder2Order(4,5) = (fu_the_u_rho - fu_the_d_rho - fd_the_u_rho + fd_the_d_rho)/(4*dthe*drho);
%Rho
numder1Order(1,5) = (fu_rho - fd_rho)/(2*drho);
numder2Order(5,5) = (fu_rho - 2*fm_ + fd_rho)/(drho^2);
numder2Order(5,1) = (fu_rho_u_nu0 - fu_rho_d_nu0 - fd_rho_u_nu0 + fd_rho_d_nu0)/(4*drho*dnu0);
numder2Order(5,2) = (fu_rho_u_kap - fu_rho_d_kap - fd_rho_u_kap + fd_rho_d_kap)/(4*drho*dkap);
numder2Order(5,3) = (fu_rho_u_eba - fu_rho_d_eba - fd_rho_u_eba + fd_rho_d_eba)/(4*drho*deba);
numder2Order(5,4) = (fu_rho_u_the - fu_rho_d_the - fd_rho_u_the + fd_rho_d_the)/(4*drho*dthe);



%%
            fprintf('******************************************VALUES and DERIVATIVES***********************************************************************************\n')
hes_e1.val; fprintf('Värde     på derivatet enligt hes_e1 : %20.12f\n',ans.val)
hes_e2.val; fprintf('Värde     på derivatet enligt hes_e2 : %20.12f\n',ans.val)
hes_e3.val; fprintf('Värde     på derivatet enligt hes_e3 : %20.12f\n',ans.val)
hes_e4.val; fprintf('Värde     på derivatet enligt hes_e4 : %20.12f\n',ans.val)
hes_e5.val; fprintf('Värde     på derivatet enligt hes_e5 : %20.12f\n',ans.val)
            fprintf('\n')
            fprintf('Derivatan på derivatet enligt grad   : '); fprintf('%20.12f  ',grad.der); fprintf('\n')
hes_e1.der; fprintf('Derivatan på derivatet enligt hes_e1 : '); fprintf('%20.12f  ',ans.val);  fprintf('\n')
hes_e2.der; fprintf('Derivatan på derivatet enligt hes_e2 : '); fprintf('%20.12f  ',ans.val);  fprintf('\n')
hes_e3.der; fprintf('Derivatan på derivatet enligt hes_e3 : '); fprintf('%20.12f  ',ans.val);  fprintf('\n')
hes_e4.der; fprintf('Derivatan på derivatet enligt hes_e4 : '); fprintf('%20.12f  ',ans.val);  fprintf('\n')
hes_e5.der; fprintf('Derivatan på derivatet enligt hes_e5 : '); fprintf('%20.12f  ',ans.val);  fprintf('\n')
            fprintf('\n')
hes_h1.der; fprintf('Derivatan på derivatet enligt hes_h1 : '); fprintf('%20.12f  ',ans.val);  fprintf('\n')
hes_h2.der; fprintf('Derivatan på derivatet enligt hes_h2 : '); fprintf('%20.12f  ',ans.val);  fprintf('\n')
hes_h3.der; fprintf('Derivatan på derivatet enligt hes_h3 : '); fprintf('%20.12f  ',ans.val);  fprintf('\n')
hes_h4.der; fprintf('Derivatan på derivatet enligt hes_h4 : '); fprintf('%20.12f  ',ans.val);  fprintf('\n')
hes_h5.der; fprintf('Derivatan på derivatet enligt hes_h5 : '); fprintf('%20.12f  ',ans.val);  fprintf('\n')
            fprintf('\n')
            fprintf('Derivatan på derivatet enligt numer  : '); fprintf('%20.12f  ',numder1Order);  fprintf('\n')
            fprintf('\n')
hes_vh.der; fprintf('2nd deri. på derivatet enligt hes_vh : '); fprintf('%20.12f  ',ans.der(1,:));  fprintf('\n')            
hes_vh.der; fprintf('                                       '); fprintf('%20.12f  ',ans.der(2,:));  fprintf('\n')            
hes_vh.der; fprintf('                                       '); fprintf('%20.12f  ',ans.der(3,:));  fprintf('\n')            
hes_vh.der; fprintf('                                       '); fprintf('%20.12f  ',ans.der(4,:));  fprintf('\n')            
hes_vh.der; fprintf('                                       '); fprintf('%20.12f  ',ans.der(5,:));  fprintf('\n')            
            fprintf('\n')
            fprintf('2nd deri. på derivatet enligt numer  : '); fprintf('%20.12f  ',numder2Order(1,:));  fprintf('\n')            
            fprintf('                                       '); fprintf('%20.12f  ',numder2Order(2,:));  fprintf('\n')            
            fprintf('                                       '); fprintf('%20.12f  ',numder2Order(3,:));  fprintf('\n')            
            fprintf('                                       '); fprintf('%20.12f  ',numder2Order(4,:));  fprintf('\n')            
            fprintf('                                       '); fprintf('%20.12f  ',numder2Order(5,:));  fprintf('\n')         
            
            fprintf('******************************************DIFFERENCES**********************************************************************************************\n')
            fprintf('Difference in first order derivateive: '); fprintf('%20.12e  ',grad.der-numder1Order);  fprintf('\n')            
            fprintf('\n')
hes_vh.der; fprintf('Difference in second order derivatvie: '); fprintf('%20.12e  ',ans.der(1,:)-numder2Order(1,:));  fprintf('\n')            
hes_vh.der; fprintf('                                       '); fprintf('%20.12e  ',ans.der(2,:)-numder2Order(2,:));  fprintf('\n')            
hes_vh.der; fprintf('                                       '); fprintf('%20.12e  ',ans.der(3,:)-numder2Order(3,:));  fprintf('\n')            
hes_vh.der; fprintf('                                       '); fprintf('%20.12e  ',ans.der(4,:)-numder2Order(4,:));  fprintf('\n')            
hes_vh.der; fprintf('                                       '); fprintf('%20.12e  ',ans.der(5,:)-numder2Order(5,:));  fprintf('\n')               
            
            
            
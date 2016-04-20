function [T,RSNP_V,RSNP_DS00iNa_m,RSNP_DS00iNa_h,RSNP_DS00iDR_n,RSNP_DS00KDyn_ko,RSNP_TW03iH_m]=solve_ode
% ------------------------------------------------------------
% Parameters:
% ------------------------------------------------------------
p=load('params.mat','p'); p=p.p;
downsample_factor=p.downsample_factor;
dt=p.dt;
T=(p.tspan(1):dt:p.tspan(2))';
ntime=length(T);
nsamp=length(1:downsample_factor:ntime);
% ------------------------------------------------------------
% Initial conditions:
% ------------------------------------------------------------
% seed the random number generator
rng(p.random_seed);
t=0; k=1;
% STATE_VARIABLES:
RSNP_V = zeros(nsamp,p.RSNP_Npop);
  RSNP_V(1,:) = -65*ones(1,p.RSNP_Npop);
RSNP_DS00iNa_m = zeros(nsamp,p.RSNP_Npop);
  RSNP_DS00iNa_m(1,:) = (((((-.2816*(( ((abs(-65-p.RSNP_DS00iNa_amV0)<p.RSNP_DS00iNa_eps).*p.RSNP_DS00iNa_eps+(abs(-65-p.RSNP_DS00iNa_amV0)>=p.RSNP_DS00iNa_eps).*-65-p.RSNP_DS00iNa_amV0))))./(-1+exp(-(( ((abs(-65-p.RSNP_DS00iNa_amV0)<p.RSNP_DS00iNa_eps).*p.RSNP_DS00iNa_eps+(abs(-65-p.RSNP_DS00iNa_amV0)>=p.RSNP_DS00iNa_eps).*-65-p.RSNP_DS00iNa_amV0)))/9.3))))./((((-.2816*(( ((abs(-65-p.RSNP_DS00iNa_amV0)<p.RSNP_DS00iNa_eps).*p.RSNP_DS00iNa_eps+(abs(-65-p.RSNP_DS00iNa_amV0)>=p.RSNP_DS00iNa_eps).*-65-p.RSNP_DS00iNa_amV0))))./(-1+exp(-(( ((abs(-65-p.RSNP_DS00iNa_amV0)<p.RSNP_DS00iNa_eps).*p.RSNP_DS00iNa_eps+(abs(-65-p.RSNP_DS00iNa_amV0)>=p.RSNP_DS00iNa_eps).*-65-p.RSNP_DS00iNa_amV0)))/9.3))))+(((.2464*(( ((abs(-65-p.RSNP_DS00iNa_bmV0)<p.RSNP_DS00iNa_eps).*p.RSNP_DS00iNa_eps+(abs(-65-p.RSNP_DS00iNa_bmV0)>=p.RSNP_DS00iNa_eps).*-65-p.RSNP_DS00iNa_bmV0))))./(-1+exp((( ((abs(-65-p.RSNP_DS00iNa_bmV0)<p.RSNP_DS00iNa_eps).*p.RSNP_DS00iNa_eps+(abs(-65-p.RSNP_DS00iNa_bmV0)>=p.RSNP_DS00iNa_eps).*-65-p.RSNP_DS00iNa_bmV0)))/6)))))))+p.RSNP_DS00iNa_IC_noise.*rand(1,p.RSNP_Npop);
RSNP_DS00iNa_h = zeros(nsamp,p.RSNP_Npop);
  RSNP_DS00iNa_h(1,:) = ((((p.RSNP_DS00iNa_hnascale.*.098./exp((-65-p.RSNP_DS00iNa_ahV0)/20)))./(((p.RSNP_DS00iNa_hnascale.*.098./exp((-65-p.RSNP_DS00iNa_ahV0)/20)))+((p.RSNP_DS00iNa_hnascale.*1.4./(1+exp(-(-65-p.RSNP_DS00iNa_bhV0)/10)))))))+p.RSNP_DS00iNa_IC_noise.*rand(1,p.RSNP_Npop);
RSNP_DS00iDR_n = zeros(nsamp,p.RSNP_Npop);
  RSNP_DS00iDR_n(1,:) = ((((-.018*(( ((abs(-65-p.RSNP_DS00iDR_anV0)<p.RSNP_DS00iDR_eps).*p.RSNP_DS00iDR_eps+(abs(-65-p.RSNP_DS00iDR_anV0)>=p.RSNP_DS00iDR_eps).*-65-p.RSNP_DS00iDR_anV0)))./(-1+exp(-(( ((abs(-65-p.RSNP_DS00iDR_anV0)<p.RSNP_DS00iDR_eps).*p.RSNP_DS00iDR_eps+(abs(-65-p.RSNP_DS00iDR_anV0)>=p.RSNP_DS00iDR_eps).*-65-p.RSNP_DS00iDR_anV0)))/25))))./(((-.018*(( ((abs(-65-p.RSNP_DS00iDR_anV0)<p.RSNP_DS00iDR_eps).*p.RSNP_DS00iDR_eps+(abs(-65-p.RSNP_DS00iDR_anV0)>=p.RSNP_DS00iDR_eps).*-65-p.RSNP_DS00iDR_anV0)))./(-1+exp(-(( ((abs(-65-p.RSNP_DS00iDR_anV0)<p.RSNP_DS00iDR_eps).*p.RSNP_DS00iDR_eps+(abs(-65-p.RSNP_DS00iDR_anV0)>=p.RSNP_DS00iDR_eps).*-65-p.RSNP_DS00iDR_anV0)))/25))))+((.0054*(( ((abs(-65-p.RSNP_DS00iDR_bnV0)<p.RSNP_DS00iDR_eps).*p.RSNP_DS00iDR_eps+(abs(-65-p.RSNP_DS00iDR_bnV0)>=p.RSNP_DS00iDR_eps).*-65-p.RSNP_DS00iDR_bnV0)))./(-1+exp((( ((abs(-65-p.RSNP_DS00iDR_bnV0)<p.RSNP_DS00iDR_eps).*p.RSNP_DS00iDR_eps+(abs(-65-p.RSNP_DS00iDR_bnV0)>=p.RSNP_DS00iDR_eps).*-65-p.RSNP_DS00iDR_bnV0)))/12)))))))+p.RSNP_DS00iDR_IC_noise.*rand(1,p.RSNP_Npop);
RSNP_DS00KDyn_ko = zeros(nsamp,p.RSNP_Npop);
  RSNP_DS00KDyn_ko(1,:) =  p.RSNP_DS00KDyn_koinf+p.RSNP_DS00KDyn_IC_noise.*rand(1,p.RSNP_Npop);
RSNP_TW03iH_m = zeros(nsamp,p.RSNP_Npop);
  RSNP_TW03iH_m(1,:) =  p.RSNP_TW03iH_IC+p.RSNP_TW03iH_IC_noise.*rand(1,p.RSNP_Npop);
% ###########################################################
% Numerical integration:
% ###########################################################
n=2;
for k=2:ntime
  t=T(k-1);
  RSNP_V_k1=(((-((p.RSNP_DS00iNa_gnaf.*RSNP_DS00iNa_m(n-1).^3.*RSNP_DS00iNa_h(n-1).*(RSNP_V(n-1)-p.RSNP_DS00iNa_ENa))))+((-((p.RSNP_DS00iDR_gkdr.*RSNP_DS00iDR_n(n-1).^4.*(RSNP_V(n-1)-((25*log(RSNP_DS00KDyn_ko(n-1)/p.RSNP_DS00iDR_ki)))))))+((-(( p.RSNP_TW03iH_gH.*RSNP_TW03iH_m(n-1).*(RSNP_V(n-1)-p.RSNP_TW03iH_E_AR))))+((-((p.RSNP_pas_gpas.*(RSNP_V(n-1)-p.RSNP_pas_epas))))))))+p.RSNP_Iapp*(t>50&t<150))./p.RSNP_Cm;
  RSNP_DS00iNa_m_k1=((((((-.2816*(( ((abs(RSNP_V(n-1)-p.RSNP_DS00iNa_amV0)<p.RSNP_DS00iNa_eps).*p.RSNP_DS00iNa_eps+(abs(RSNP_V(n-1)-p.RSNP_DS00iNa_amV0)>=p.RSNP_DS00iNa_eps).*RSNP_V(n-1)-p.RSNP_DS00iNa_amV0))))./(-1+exp(-(( ((abs(RSNP_V(n-1)-p.RSNP_DS00iNa_amV0)<p.RSNP_DS00iNa_eps).*p.RSNP_DS00iNa_eps+(abs(RSNP_V(n-1)-p.RSNP_DS00iNa_amV0)>=p.RSNP_DS00iNa_eps).*RSNP_V(n-1)-p.RSNP_DS00iNa_amV0)))/9.3))))./((((-.2816*(( ((abs(RSNP_V(n-1)-p.RSNP_DS00iNa_amV0)<p.RSNP_DS00iNa_eps).*p.RSNP_DS00iNa_eps+(abs(RSNP_V(n-1)-p.RSNP_DS00iNa_amV0)>=p.RSNP_DS00iNa_eps).*RSNP_V(n-1)-p.RSNP_DS00iNa_amV0))))./(-1+exp(-(( ((abs(RSNP_V(n-1)-p.RSNP_DS00iNa_amV0)<p.RSNP_DS00iNa_eps).*p.RSNP_DS00iNa_eps+(abs(RSNP_V(n-1)-p.RSNP_DS00iNa_amV0)>=p.RSNP_DS00iNa_eps).*RSNP_V(n-1)-p.RSNP_DS00iNa_amV0)))/9.3))))+(((.2464*(( ((abs(RSNP_V(n-1)-p.RSNP_DS00iNa_bmV0)<p.RSNP_DS00iNa_eps).*p.RSNP_DS00iNa_eps+(abs(RSNP_V(n-1)-p.RSNP_DS00iNa_bmV0)>=p.RSNP_DS00iNa_eps).*RSNP_V(n-1)-p.RSNP_DS00iNa_bmV0))))./(-1+exp((( ((abs(RSNP_V(n-1)-p.RSNP_DS00iNa_bmV0)<p.RSNP_DS00iNa_eps).*p.RSNP_DS00iNa_eps+(abs(RSNP_V(n-1)-p.RSNP_DS00iNa_bmV0)>=p.RSNP_DS00iNa_eps).*RSNP_V(n-1)-p.RSNP_DS00iNa_bmV0)))/6)))))))-RSNP_DS00iNa_m(n-1))./((1./((((-.2816*(( ((abs(RSNP_V(n-1)-p.RSNP_DS00iNa_amV0)<p.RSNP_DS00iNa_eps).*p.RSNP_DS00iNa_eps+(abs(RSNP_V(n-1)-p.RSNP_DS00iNa_amV0)>=p.RSNP_DS00iNa_eps).*RSNP_V(n-1)-p.RSNP_DS00iNa_amV0))))./(-1+exp(-(( ((abs(RSNP_V(n-1)-p.RSNP_DS00iNa_amV0)<p.RSNP_DS00iNa_eps).*p.RSNP_DS00iNa_eps+(abs(RSNP_V(n-1)-p.RSNP_DS00iNa_amV0)>=p.RSNP_DS00iNa_eps).*RSNP_V(n-1)-p.RSNP_DS00iNa_amV0)))/9.3))))+(((.2464*(( ((abs(RSNP_V(n-1)-p.RSNP_DS00iNa_bmV0)<p.RSNP_DS00iNa_eps).*p.RSNP_DS00iNa_eps+(abs(RSNP_V(n-1)-p.RSNP_DS00iNa_bmV0)>=p.RSNP_DS00iNa_eps).*RSNP_V(n-1)-p.RSNP_DS00iNa_bmV0))))./(-1+exp((( ((abs(RSNP_V(n-1)-p.RSNP_DS00iNa_bmV0)<p.RSNP_DS00iNa_eps).*p.RSNP_DS00iNa_eps+(abs(RSNP_V(n-1)-p.RSNP_DS00iNa_bmV0)>=p.RSNP_DS00iNa_eps).*RSNP_V(n-1)-p.RSNP_DS00iNa_bmV0)))/6)))))));
  RSNP_DS00iNa_h_k1=(((((p.RSNP_DS00iNa_hnascale.*.098./exp((RSNP_V(n-1)-p.RSNP_DS00iNa_ahV0)/20)))./(((p.RSNP_DS00iNa_hnascale.*.098./exp((RSNP_V(n-1)-p.RSNP_DS00iNa_ahV0)/20)))+((p.RSNP_DS00iNa_hnascale.*1.4./(1+exp(-(RSNP_V(n-1)-p.RSNP_DS00iNa_bhV0)/10)))))))-RSNP_DS00iNa_h(n-1))./((1./(((p.RSNP_DS00iNa_hnascale.*.098./exp((RSNP_V(n-1)-p.RSNP_DS00iNa_ahV0)/20)))+((p.RSNP_DS00iNa_hnascale.*1.4./(1+exp(-(RSNP_V(n-1)-p.RSNP_DS00iNa_bhV0)/10)))))));
  RSNP_DS00iDR_n_k1=(((((-.018*(( ((abs(RSNP_V(n-1)-p.RSNP_DS00iDR_anV0)<p.RSNP_DS00iDR_eps).*p.RSNP_DS00iDR_eps+(abs(RSNP_V(n-1)-p.RSNP_DS00iDR_anV0)>=p.RSNP_DS00iDR_eps).*RSNP_V(n-1)-p.RSNP_DS00iDR_anV0)))./(-1+exp(-(( ((abs(RSNP_V(n-1)-p.RSNP_DS00iDR_anV0)<p.RSNP_DS00iDR_eps).*p.RSNP_DS00iDR_eps+(abs(RSNP_V(n-1)-p.RSNP_DS00iDR_anV0)>=p.RSNP_DS00iDR_eps).*RSNP_V(n-1)-p.RSNP_DS00iDR_anV0)))/25))))./(((-.018*(( ((abs(RSNP_V(n-1)-p.RSNP_DS00iDR_anV0)<p.RSNP_DS00iDR_eps).*p.RSNP_DS00iDR_eps+(abs(RSNP_V(n-1)-p.RSNP_DS00iDR_anV0)>=p.RSNP_DS00iDR_eps).*RSNP_V(n-1)-p.RSNP_DS00iDR_anV0)))./(-1+exp(-(( ((abs(RSNP_V(n-1)-p.RSNP_DS00iDR_anV0)<p.RSNP_DS00iDR_eps).*p.RSNP_DS00iDR_eps+(abs(RSNP_V(n-1)-p.RSNP_DS00iDR_anV0)>=p.RSNP_DS00iDR_eps).*RSNP_V(n-1)-p.RSNP_DS00iDR_anV0)))/25))))+((.0054*(( ((abs(RSNP_V(n-1)-p.RSNP_DS00iDR_bnV0)<p.RSNP_DS00iDR_eps).*p.RSNP_DS00iDR_eps+(abs(RSNP_V(n-1)-p.RSNP_DS00iDR_bnV0)>=p.RSNP_DS00iDR_eps).*RSNP_V(n-1)-p.RSNP_DS00iDR_bnV0)))./(-1+exp((( ((abs(RSNP_V(n-1)-p.RSNP_DS00iDR_bnV0)<p.RSNP_DS00iDR_eps).*p.RSNP_DS00iDR_eps+(abs(RSNP_V(n-1)-p.RSNP_DS00iDR_bnV0)>=p.RSNP_DS00iDR_eps).*RSNP_V(n-1)-p.RSNP_DS00iDR_bnV0)))/12)))))))-RSNP_DS00iDR_n(n-1))./((1./(((-.018*(( ((abs(RSNP_V(n-1)-p.RSNP_DS00iDR_anV0)<p.RSNP_DS00iDR_eps).*p.RSNP_DS00iDR_eps+(abs(RSNP_V(n-1)-p.RSNP_DS00iDR_anV0)>=p.RSNP_DS00iDR_eps).*RSNP_V(n-1)-p.RSNP_DS00iDR_anV0)))./(-1+exp(-(( ((abs(RSNP_V(n-1)-p.RSNP_DS00iDR_anV0)<p.RSNP_DS00iDR_eps).*p.RSNP_DS00iDR_eps+(abs(RSNP_V(n-1)-p.RSNP_DS00iDR_anV0)>=p.RSNP_DS00iDR_eps).*RSNP_V(n-1)-p.RSNP_DS00iDR_anV0)))/25))))+((.0054*(( ((abs(RSNP_V(n-1)-p.RSNP_DS00iDR_bnV0)<p.RSNP_DS00iDR_eps).*p.RSNP_DS00iDR_eps+(abs(RSNP_V(n-1)-p.RSNP_DS00iDR_bnV0)>=p.RSNP_DS00iDR_eps).*RSNP_V(n-1)-p.RSNP_DS00iDR_bnV0)))./(-1+exp((( ((abs(RSNP_V(n-1)-p.RSNP_DS00iDR_bnV0)<p.RSNP_DS00iDR_eps).*p.RSNP_DS00iDR_eps+(abs(RSNP_V(n-1)-p.RSNP_DS00iDR_bnV0)>=p.RSNP_DS00iDR_eps).*RSNP_V(n-1)-p.RSNP_DS00iDR_bnV0)))/12)))))));
  RSNP_DS00KDyn_ko_k1= p.RSNP_DS00KDyn_KAF.*((((p.RSNP_DS00iDR_gkdr.*RSNP_DS00iDR_n(n-1).^4.*(RSNP_V(n-1)-((25*log(RSNP_DS00KDyn_ko(n-1)/p.RSNP_DS00iDR_ki))))))))./(p.RSNP_DS00KDyn_faraday*p.RSNP_DS00KDyn_VshellK)+(p.RSNP_DS00KDyn_koinf-RSNP_DS00KDyn_ko(n-1))./p.RSNP_DS00KDyn_tauK;
  RSNP_TW03iH_m_k1= (( p.RSNP_TW03iH_c_ARaM.*((( 1 ./ (1+exp((p.RSNP_TW03iH_AR_V12-RSNP_V(n-1))/p.RSNP_TW03iH_AR_k)))) ./ (( 1./(p.RSNP_TW03iH_AR_L.*exp(-14.6-.086*RSNP_V(n-1))+p.RSNP_TW03iH_AR_R.*exp(-1.87+.07*RSNP_V(n-1)))))))).*(1-RSNP_TW03iH_m(n-1))-(( p.RSNP_TW03iH_c_ARbM.*((1-(( 1 ./ (1+exp((p.RSNP_TW03iH_AR_V12-RSNP_V(n-1))/p.RSNP_TW03iH_AR_k)))))./(( 1./(p.RSNP_TW03iH_AR_L.*exp(-14.6-.086*RSNP_V(n-1))+p.RSNP_TW03iH_AR_R.*exp(-1.87+.07*RSNP_V(n-1)))))))).*RSNP_TW03iH_m(n-1);
  % ------------------------------------------------------------
  % Update state variables:
  % ------------------------------------------------------------
  RSNP_V(n)=RSNP_V(n-1)+dt*RSNP_V_k1;
  RSNP_DS00iNa_m(n)=RSNP_DS00iNa_m(n-1)+dt*RSNP_DS00iNa_m_k1;
  RSNP_DS00iNa_h(n)=RSNP_DS00iNa_h(n-1)+dt*RSNP_DS00iNa_h_k1;
  RSNP_DS00iDR_n(n)=RSNP_DS00iDR_n(n-1)+dt*RSNP_DS00iDR_n_k1;
  RSNP_DS00KDyn_ko(n)=RSNP_DS00KDyn_ko(n-1)+dt*RSNP_DS00KDyn_ko_k1;
  RSNP_TW03iH_m(n)=RSNP_TW03iH_m(n-1)+dt*RSNP_TW03iH_m_k1;
  n=n+1;
end
T=T(1:downsample_factor:ntime);

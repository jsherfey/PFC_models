function [T,FS_V,FS_DS00iNa_m,FS_DS00iNa_h,FS_DS00iDR_n,FS_DS00KDyn_ko]=solve_ode
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
FS_V = zeros(nsamp,p.FS_Npop);
  FS_V(1,:) = -65*ones(1,p.FS_Npop);
FS_DS00iNa_m = zeros(nsamp,p.FS_Npop);
  FS_DS00iNa_m(1,:) = (((((-.2816*(( ((abs(-65-p.FS_DS00iNa_amV0)<p.FS_DS00iNa_eps).*p.FS_DS00iNa_eps+(abs(-65-p.FS_DS00iNa_amV0)>=p.FS_DS00iNa_eps).*-65-p.FS_DS00iNa_amV0))))./(-1+exp(-(( ((abs(-65-p.FS_DS00iNa_amV0)<p.FS_DS00iNa_eps).*p.FS_DS00iNa_eps+(abs(-65-p.FS_DS00iNa_amV0)>=p.FS_DS00iNa_eps).*-65-p.FS_DS00iNa_amV0)))/9.3))))./((((-.2816*(( ((abs(-65-p.FS_DS00iNa_amV0)<p.FS_DS00iNa_eps).*p.FS_DS00iNa_eps+(abs(-65-p.FS_DS00iNa_amV0)>=p.FS_DS00iNa_eps).*-65-p.FS_DS00iNa_amV0))))./(-1+exp(-(( ((abs(-65-p.FS_DS00iNa_amV0)<p.FS_DS00iNa_eps).*p.FS_DS00iNa_eps+(abs(-65-p.FS_DS00iNa_amV0)>=p.FS_DS00iNa_eps).*-65-p.FS_DS00iNa_amV0)))/9.3))))+(((.2464*(( ((abs(-65-p.FS_DS00iNa_bmV0)<p.FS_DS00iNa_eps).*p.FS_DS00iNa_eps+(abs(-65-p.FS_DS00iNa_bmV0)>=p.FS_DS00iNa_eps).*-65-p.FS_DS00iNa_bmV0))))./(-1+exp((( ((abs(-65-p.FS_DS00iNa_bmV0)<p.FS_DS00iNa_eps).*p.FS_DS00iNa_eps+(abs(-65-p.FS_DS00iNa_bmV0)>=p.FS_DS00iNa_eps).*-65-p.FS_DS00iNa_bmV0)))/6)))))))+p.FS_DS00iNa_IC_noise.*rand(1,p.FS_Npop);
FS_DS00iNa_h = zeros(nsamp,p.FS_Npop);
  FS_DS00iNa_h(1,:) = ((((p.FS_DS00iNa_hnascale.*.098./exp((-65-p.FS_DS00iNa_ahV0)/20)))./(((p.FS_DS00iNa_hnascale.*.098./exp((-65-p.FS_DS00iNa_ahV0)/20)))+((p.FS_DS00iNa_hnascale.*1.4./(1+exp(-(-65-p.FS_DS00iNa_bhV0)/10)))))))+p.FS_DS00iNa_IC_noise.*rand(1,p.FS_Npop);
FS_DS00iDR_n = zeros(nsamp,p.FS_Npop);
  FS_DS00iDR_n(1,:) = ((((-.018*(( ((abs(-65-p.FS_DS00iDR_anV0)<p.FS_DS00iDR_eps).*p.FS_DS00iDR_eps+(abs(-65-p.FS_DS00iDR_anV0)>=p.FS_DS00iDR_eps).*-65-p.FS_DS00iDR_anV0)))./(-1+exp(-(( ((abs(-65-p.FS_DS00iDR_anV0)<p.FS_DS00iDR_eps).*p.FS_DS00iDR_eps+(abs(-65-p.FS_DS00iDR_anV0)>=p.FS_DS00iDR_eps).*-65-p.FS_DS00iDR_anV0)))/25))))./(((-.018*(( ((abs(-65-p.FS_DS00iDR_anV0)<p.FS_DS00iDR_eps).*p.FS_DS00iDR_eps+(abs(-65-p.FS_DS00iDR_anV0)>=p.FS_DS00iDR_eps).*-65-p.FS_DS00iDR_anV0)))./(-1+exp(-(( ((abs(-65-p.FS_DS00iDR_anV0)<p.FS_DS00iDR_eps).*p.FS_DS00iDR_eps+(abs(-65-p.FS_DS00iDR_anV0)>=p.FS_DS00iDR_eps).*-65-p.FS_DS00iDR_anV0)))/25))))+((.0054*(( ((abs(-65-p.FS_DS00iDR_bnV0)<p.FS_DS00iDR_eps).*p.FS_DS00iDR_eps+(abs(-65-p.FS_DS00iDR_bnV0)>=p.FS_DS00iDR_eps).*-65-p.FS_DS00iDR_bnV0)))./(-1+exp((( ((abs(-65-p.FS_DS00iDR_bnV0)<p.FS_DS00iDR_eps).*p.FS_DS00iDR_eps+(abs(-65-p.FS_DS00iDR_bnV0)>=p.FS_DS00iDR_eps).*-65-p.FS_DS00iDR_bnV0)))/12)))))))+p.FS_DS00iDR_IC_noise.*rand(1,p.FS_Npop);
FS_DS00KDyn_ko = zeros(nsamp,p.FS_Npop);
  FS_DS00KDyn_ko(1,:) =  p.FS_DS00KDyn_koinf+p.FS_DS00KDyn_IC_noise.*rand(1,p.FS_Npop);
% ###########################################################
% Numerical integration:
% ###########################################################
n=2;
for k=2:ntime
  t=T(k-1);
  FS_V_k1=(((-((p.FS_DS00iNa_gnaf.*FS_DS00iNa_m(n-1).^3.*FS_DS00iNa_h(n-1).*(FS_V(n-1)-p.FS_DS00iNa_ENa))))+((-((p.FS_DS00iDR_gkdr.*FS_DS00iDR_n(n-1).^4.*(FS_V(n-1)-((25*log(FS_DS00KDyn_ko(n-1)/p.FS_DS00iDR_ki)))))))+((-((p.FS_pas_gpas.*(FS_V(n-1)-p.FS_pas_epas)))))))+p.FS_Iapp)./p.FS_Cm;
  FS_DS00iNa_m_k1=((((((-.2816*(( ((abs(FS_V(n-1)-p.FS_DS00iNa_amV0)<p.FS_DS00iNa_eps).*p.FS_DS00iNa_eps+(abs(FS_V(n-1)-p.FS_DS00iNa_amV0)>=p.FS_DS00iNa_eps).*FS_V(n-1)-p.FS_DS00iNa_amV0))))./(-1+exp(-(( ((abs(FS_V(n-1)-p.FS_DS00iNa_amV0)<p.FS_DS00iNa_eps).*p.FS_DS00iNa_eps+(abs(FS_V(n-1)-p.FS_DS00iNa_amV0)>=p.FS_DS00iNa_eps).*FS_V(n-1)-p.FS_DS00iNa_amV0)))/9.3))))./((((-.2816*(( ((abs(FS_V(n-1)-p.FS_DS00iNa_amV0)<p.FS_DS00iNa_eps).*p.FS_DS00iNa_eps+(abs(FS_V(n-1)-p.FS_DS00iNa_amV0)>=p.FS_DS00iNa_eps).*FS_V(n-1)-p.FS_DS00iNa_amV0))))./(-1+exp(-(( ((abs(FS_V(n-1)-p.FS_DS00iNa_amV0)<p.FS_DS00iNa_eps).*p.FS_DS00iNa_eps+(abs(FS_V(n-1)-p.FS_DS00iNa_amV0)>=p.FS_DS00iNa_eps).*FS_V(n-1)-p.FS_DS00iNa_amV0)))/9.3))))+(((.2464*(( ((abs(FS_V(n-1)-p.FS_DS00iNa_bmV0)<p.FS_DS00iNa_eps).*p.FS_DS00iNa_eps+(abs(FS_V(n-1)-p.FS_DS00iNa_bmV0)>=p.FS_DS00iNa_eps).*FS_V(n-1)-p.FS_DS00iNa_bmV0))))./(-1+exp((( ((abs(FS_V(n-1)-p.FS_DS00iNa_bmV0)<p.FS_DS00iNa_eps).*p.FS_DS00iNa_eps+(abs(FS_V(n-1)-p.FS_DS00iNa_bmV0)>=p.FS_DS00iNa_eps).*FS_V(n-1)-p.FS_DS00iNa_bmV0)))/6)))))))-FS_DS00iNa_m(n-1))./((1./((((-.2816*(( ((abs(FS_V(n-1)-p.FS_DS00iNa_amV0)<p.FS_DS00iNa_eps).*p.FS_DS00iNa_eps+(abs(FS_V(n-1)-p.FS_DS00iNa_amV0)>=p.FS_DS00iNa_eps).*FS_V(n-1)-p.FS_DS00iNa_amV0))))./(-1+exp(-(( ((abs(FS_V(n-1)-p.FS_DS00iNa_amV0)<p.FS_DS00iNa_eps).*p.FS_DS00iNa_eps+(abs(FS_V(n-1)-p.FS_DS00iNa_amV0)>=p.FS_DS00iNa_eps).*FS_V(n-1)-p.FS_DS00iNa_amV0)))/9.3))))+(((.2464*(( ((abs(FS_V(n-1)-p.FS_DS00iNa_bmV0)<p.FS_DS00iNa_eps).*p.FS_DS00iNa_eps+(abs(FS_V(n-1)-p.FS_DS00iNa_bmV0)>=p.FS_DS00iNa_eps).*FS_V(n-1)-p.FS_DS00iNa_bmV0))))./(-1+exp((( ((abs(FS_V(n-1)-p.FS_DS00iNa_bmV0)<p.FS_DS00iNa_eps).*p.FS_DS00iNa_eps+(abs(FS_V(n-1)-p.FS_DS00iNa_bmV0)>=p.FS_DS00iNa_eps).*FS_V(n-1)-p.FS_DS00iNa_bmV0)))/6)))))));
  FS_DS00iNa_h_k1=(((((p.FS_DS00iNa_hnascale.*.098./exp((FS_V(n-1)-p.FS_DS00iNa_ahV0)/20)))./(((p.FS_DS00iNa_hnascale.*.098./exp((FS_V(n-1)-p.FS_DS00iNa_ahV0)/20)))+((p.FS_DS00iNa_hnascale.*1.4./(1+exp(-(FS_V(n-1)-p.FS_DS00iNa_bhV0)/10)))))))-FS_DS00iNa_h(n-1))./((1./(((p.FS_DS00iNa_hnascale.*.098./exp((FS_V(n-1)-p.FS_DS00iNa_ahV0)/20)))+((p.FS_DS00iNa_hnascale.*1.4./(1+exp(-(FS_V(n-1)-p.FS_DS00iNa_bhV0)/10)))))));
  FS_DS00iDR_n_k1=(((((-.018*(( ((abs(FS_V(n-1)-p.FS_DS00iDR_anV0)<p.FS_DS00iDR_eps).*p.FS_DS00iDR_eps+(abs(FS_V(n-1)-p.FS_DS00iDR_anV0)>=p.FS_DS00iDR_eps).*FS_V(n-1)-p.FS_DS00iDR_anV0)))./(-1+exp(-(( ((abs(FS_V(n-1)-p.FS_DS00iDR_anV0)<p.FS_DS00iDR_eps).*p.FS_DS00iDR_eps+(abs(FS_V(n-1)-p.FS_DS00iDR_anV0)>=p.FS_DS00iDR_eps).*FS_V(n-1)-p.FS_DS00iDR_anV0)))/25))))./(((-.018*(( ((abs(FS_V(n-1)-p.FS_DS00iDR_anV0)<p.FS_DS00iDR_eps).*p.FS_DS00iDR_eps+(abs(FS_V(n-1)-p.FS_DS00iDR_anV0)>=p.FS_DS00iDR_eps).*FS_V(n-1)-p.FS_DS00iDR_anV0)))./(-1+exp(-(( ((abs(FS_V(n-1)-p.FS_DS00iDR_anV0)<p.FS_DS00iDR_eps).*p.FS_DS00iDR_eps+(abs(FS_V(n-1)-p.FS_DS00iDR_anV0)>=p.FS_DS00iDR_eps).*FS_V(n-1)-p.FS_DS00iDR_anV0)))/25))))+((.0054*(( ((abs(FS_V(n-1)-p.FS_DS00iDR_bnV0)<p.FS_DS00iDR_eps).*p.FS_DS00iDR_eps+(abs(FS_V(n-1)-p.FS_DS00iDR_bnV0)>=p.FS_DS00iDR_eps).*FS_V(n-1)-p.FS_DS00iDR_bnV0)))./(-1+exp((( ((abs(FS_V(n-1)-p.FS_DS00iDR_bnV0)<p.FS_DS00iDR_eps).*p.FS_DS00iDR_eps+(abs(FS_V(n-1)-p.FS_DS00iDR_bnV0)>=p.FS_DS00iDR_eps).*FS_V(n-1)-p.FS_DS00iDR_bnV0)))/12)))))))-FS_DS00iDR_n(n-1))./((1./(((-.018*(( ((abs(FS_V(n-1)-p.FS_DS00iDR_anV0)<p.FS_DS00iDR_eps).*p.FS_DS00iDR_eps+(abs(FS_V(n-1)-p.FS_DS00iDR_anV0)>=p.FS_DS00iDR_eps).*FS_V(n-1)-p.FS_DS00iDR_anV0)))./(-1+exp(-(( ((abs(FS_V(n-1)-p.FS_DS00iDR_anV0)<p.FS_DS00iDR_eps).*p.FS_DS00iDR_eps+(abs(FS_V(n-1)-p.FS_DS00iDR_anV0)>=p.FS_DS00iDR_eps).*FS_V(n-1)-p.FS_DS00iDR_anV0)))/25))))+((.0054*(( ((abs(FS_V(n-1)-p.FS_DS00iDR_bnV0)<p.FS_DS00iDR_eps).*p.FS_DS00iDR_eps+(abs(FS_V(n-1)-p.FS_DS00iDR_bnV0)>=p.FS_DS00iDR_eps).*FS_V(n-1)-p.FS_DS00iDR_bnV0)))./(-1+exp((( ((abs(FS_V(n-1)-p.FS_DS00iDR_bnV0)<p.FS_DS00iDR_eps).*p.FS_DS00iDR_eps+(abs(FS_V(n-1)-p.FS_DS00iDR_bnV0)>=p.FS_DS00iDR_eps).*FS_V(n-1)-p.FS_DS00iDR_bnV0)))/12)))))));
  FS_DS00KDyn_ko_k1= p.FS_DS00KDyn_KAF.*((((p.FS_DS00iDR_gkdr.*FS_DS00iDR_n(n-1).^4.*(FS_V(n-1)-((25*log(FS_DS00KDyn_ko(n-1)/p.FS_DS00iDR_ki))))))))./(p.FS_DS00KDyn_faraday*p.FS_DS00KDyn_VshellK)+(p.FS_DS00KDyn_koinf-FS_DS00KDyn_ko(n-1))./p.FS_DS00KDyn_tauK;
  % ------------------------------------------------------------
  % Update state variables:
  % ------------------------------------------------------------
  FS_V(n)=FS_V(n-1)+dt*FS_V_k1;
  FS_DS00iNa_m(n)=FS_DS00iNa_m(n-1)+dt*FS_DS00iNa_m_k1;
  FS_DS00iNa_h(n)=FS_DS00iNa_h(n-1)+dt*FS_DS00iNa_h_k1;
  FS_DS00iDR_n(n)=FS_DS00iDR_n(n-1)+dt*FS_DS00iDR_n_k1;
  FS_DS00KDyn_ko(n)=FS_DS00KDyn_ko(n-1)+dt*FS_DS00KDyn_ko_k1;
  n=n+1;
end
T=T(1:downsample_factor:ntime);

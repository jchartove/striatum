function [T,soma_V,soma_somaGolombK_a,soma_somaGolombK_b,soma_somaGolombKdr_n,soma_somaGolombNa_h,dend_V,dend_dendGolombK_a,dend_dendGolombK_b,dend_dendGolombKdr_n,dend_dendGolombNa_h,D1_V,D1_naCurrentMSN_m,D1_naCurrentMSN_h,D1_kCurrentMSN_m,D1_mCurrentMSN_m,D2_V,D2_naCurrentMSN_m,D2_naCurrentMSN_h,D2_kCurrentMSN_m,D2_mCurrentMSN_m,soma_soma_somaSomaiSYN_s,D1_soma_somaMSNiSYN_s,D2_soma_somaMSNiSYN_s,D1_D1_gabaRecInputMSN_s,D2_D2_gabaRecInputMSN_s,dend_dendGolombK_gd_dend,dend_dendGolombKdr_gkdr_dend,dend_dendGolombNa_gna_dend,dend_dendInput_tonic2,dend_dendiMultiPoissonExp_Ge,dend_dendiMultiPoissonExp_Gi,D1_naCurrentMSN_V_IC,D1_kCurrentMSN_V_IC,D1_mCurrentMSN_Qs,D1_mCurrentMSN_V_IC,D1_injectedCurrentD1_freq,D2_naCurrentMSN_V_IC,D2_kCurrentMSN_V_IC,D2_mCurrentMSN_Qs,D2_mCurrentMSN_V_IC,D2_injectedCurrentD2_freq,soma_soma_somaSomaiSYN_gsyn2,soma_soma_somaSomaiSYN_mask,dend_dend_dendDendiGAP_g_GAP2,dend_dend_dendDendiGAP_mask,D1_soma_somaMSNiSYN_indegree,D1_soma_somaMSNiSYN_gsyn2,D1_soma_somaMSNiSYN_mask,D2_soma_somaMSNiSYN_indegree,D2_soma_somaMSNiSYN_gsyn2,D2_soma_somaMSNiSYN_mask,D1_D1_gabaRecInputMSN_netcon,D2_D2_gabaRecInputMSN_netcon]=solve_ode
% ------------------------------------------------------------
% Parameters:
% ------------------------------------------------------------
params = load('params.mat','p');
p = params.p;
downsample_factor=p.downsample_factor;
dt=p.dt;
T=(p.tspan(1):dt:p.tspan(2))';
ntime=length(T);
nsamp=length(1:downsample_factor:ntime);
% ------------------------------------------------------------
% Fixed variables:
% ------------------------------------------------------------
dend_dendGolombK_gd_dend =  p.dend_dendGolombK_gd/10;
dend_dendGolombKdr_gkdr_dend =  p.dend_dendGolombKdr_gkdr/10;
dend_dendGolombNa_gna_dend =  p.dend_dendGolombNa_gna/10;
dend_dendInput_tonic2 =  p.dend_dendInput_tonic + 10.*p.dend_dendInput_DA;
dend_dendiMultiPoissonExp_Ge =  corrPoisson(p.dend_Npop, p.dend_dendiMultiPoissonExp_N_einputs, p.dend_dendiMultiPoissonExp_rate, p.dend_dendiMultiPoissonExp_tau_i, p.dend_dendiMultiPoissonExp_tau_1, 2, .5, p.dend_dendiMultiPoissonExp_Tfinal, dt, p.dend_dendiMultiPoissonExp_fraction_shared);
dend_dendiMultiPoissonExp_Gi =  corrPoisson(p.dend_Npop, p.dend_dendiMultiPoissonExp_N_iinputs, p.dend_dendiMultiPoissonExp_rate, p.dend_dendiMultiPoissonExp_tau_i, p.dend_dendiMultiPoissonExp_tau_1, 5, .5, p.dend_dendiMultiPoissonExp_Tfinal, dt, p.dend_dendiMultiPoissonExp_fraction_shared);
D1_naCurrentMSN_V_IC = -63.00000000;
D1_kCurrentMSN_V_IC = -63.00000000;
D1_mCurrentMSN_Qs =  p.D1_mCurrentMSN_Q10^(.1*(37-23));
D1_mCurrentMSN_V_IC = -63.00000000;
D1_injectedCurrentD1_freq =  1/(2*pi);
D2_naCurrentMSN_V_IC = -63.00000000;
D2_kCurrentMSN_V_IC = -63.00000000;
D2_mCurrentMSN_Qs =  p.D2_mCurrentMSN_Q10^(.1*(37-23));
D2_mCurrentMSN_V_IC = -63.00000000;
D2_injectedCurrentD2_freq =  1/(2*pi);
soma_soma_somaSomaiSYN_gsyn2 =  p.soma_soma_somaSomaiSYN_gsyn - p.soma_soma_somaSomaiSYN_DA.*(p.soma_soma_somaSomaiSYN_gsyn - 0.005);
soma_soma_somaSomaiSYN_mask =  genmask(p.soma_Npop,p.soma_Npop,p.soma_soma_somaSomaiSYN_i_con,soma_soma_somaSomaiSYN_gsyn2,1,1,0);
dend_dend_dendDendiGAP_g_GAP2 =  p.dend_dend_dendDendiGAP_g_GAP + p.dend_dend_dendDendiGAP_DA.*(0.4-p.dend_dend_dendDendiGAP_g_GAP);
dend_dend_dendDendiGAP_mask =  genmask(p.dend_Npop,p.dend_Npop,p.dend_dend_dendDendiGAP_gcon,dend_dend_dendDendiGAP_g_GAP2,0,1,0);
D1_soma_somaMSNiSYN_indegree =  p.D1_soma_somaMSNiSYN_i_con*(p.soma_Npop-p.D1_soma_somaMSNiSYN_ko);
D1_soma_somaMSNiSYN_gsyn2 =  (10*p.D1_soma_somaMSNiSYN_gsyn)/(D1_soma_somaMSNiSYN_indegree+1);
D1_soma_somaMSNiSYN_mask =  genmask(p.soma_Npop,p.D1_Npop,p.D1_soma_somaMSNiSYN_i_con,D1_soma_somaMSNiSYN_gsyn2,1,0,p.D1_soma_somaMSNiSYN_ko);
D2_soma_somaMSNiSYN_indegree =  p.D2_soma_somaMSNiSYN_i_con*(p.soma_Npop-p.D2_soma_somaMSNiSYN_ko);
D2_soma_somaMSNiSYN_gsyn2 =  (10*p.D2_soma_somaMSNiSYN_gsyn)/(D2_soma_somaMSNiSYN_indegree+1);
D2_soma_somaMSNiSYN_mask =  genmask(p.soma_Npop,p.D2_Npop,p.D2_soma_somaMSNiSYN_i_con,D2_soma_somaMSNiSYN_gsyn2,1,0,p.D2_soma_somaMSNiSYN_ko);
D1_D1_gabaRecInputMSN_netcon =  ones(p.D1_Npop)-eye(p.D1_Npop);
D2_D2_gabaRecInputMSN_netcon =  ones(p.D2_Npop)-eye(p.D2_Npop);
% ------------------------------------------------------------
% Initial conditions:
% ------------------------------------------------------------
% seed the random number generator
rng(p.random_seed);
t=0; k=1;

% STATE_VARIABLES:
soma_V = zeros(nsamp,p.soma_Npop);
  soma_V(1,:) = zeros(1,p.soma_Npop);
soma_somaGolombK_a = zeros(nsamp,p.soma_Npop);
  soma_somaGolombK_a(1,:) =  0.15 + 0.60*rand(1,p.soma_Npop);
soma_somaGolombK_b = zeros(nsamp,p.soma_Npop);
  soma_somaGolombK_b(1,:) =  0.15 + 0.45*rand(1,p.soma_Npop);
soma_somaGolombKdr_n = zeros(nsamp,p.soma_Npop);
  soma_somaGolombKdr_n(1,:) =  0.05 + 0.45*rand(1,p.soma_Npop);
soma_somaGolombNa_h = zeros(nsamp,p.soma_Npop);
  soma_somaGolombNa_h(1,:) =  0.05 + 0.85*rand(1,p.soma_Npop);
dend_V = zeros(nsamp,p.dend_Npop);
  dend_V(1,:) = zeros(1,p.dend_Npop);
dend_dendGolombK_a = zeros(nsamp,p.dend_Npop);
  dend_dendGolombK_a(1,:) =  0.15 + 0.60*rand(1,p.dend_Npop);
dend_dendGolombK_b = zeros(nsamp,p.dend_Npop);
  dend_dendGolombK_b(1,:) =  0.15 + 0.45*rand(1,p.dend_Npop);
dend_dendGolombKdr_n = zeros(nsamp,p.dend_Npop);
  dend_dendGolombKdr_n(1,:) =  0.05 + 0.45*rand(1,p.dend_Npop);
dend_dendGolombNa_h = zeros(nsamp,p.dend_Npop);
  dend_dendGolombNa_h(1,:) =  0.05 + 0.85*rand(1,p.dend_Npop);
D1_V = zeros(nsamp,p.D1_Npop);
  D1_V(1,:) = zeros(1,p.D1_Npop);
D1_naCurrentMSN_m = zeros(nsamp,p.D1_Npop);
  D1_naCurrentMSN_m(1,:) =  0.32*(D1_naCurrentMSN_V_IC+54)/(1-exp(-(D1_naCurrentMSN_V_IC+54)/4))/(0.32*(D1_naCurrentMSN_V_IC+54)/(1-exp(-(D1_naCurrentMSN_V_IC+54)/4))+0.28*(D1_naCurrentMSN_V_IC+27)/(exp((D1_naCurrentMSN_V_IC+27)/5)-1))*randn(1,p.D1_Npop);
D1_naCurrentMSN_h = zeros(nsamp,p.D1_Npop);
  D1_naCurrentMSN_h(1,:) =  0.128*exp(-(D1_naCurrentMSN_V_IC+50)/18)/(0.128*exp(-(D1_naCurrentMSN_V_IC+50)/18)+4/(1+exp(-(D1_naCurrentMSN_V_IC+27)/5)))*randn(1,p.D1_Npop);
D1_kCurrentMSN_m = zeros(nsamp,p.D1_Npop);
  D1_kCurrentMSN_m(1,:) =  0.032*(D1_kCurrentMSN_V_IC+52)/(1-exp(-(D1_kCurrentMSN_V_IC+52)/5))/(0.032*(D1_kCurrentMSN_V_IC+52)/(1-exp(-(D1_kCurrentMSN_V_IC+52)/5))+0.5*exp(-(D1_kCurrentMSN_V_IC+57)/40))*randn(1,p.D1_Npop);
D1_mCurrentMSN_m = zeros(nsamp,p.D1_Npop);
  D1_mCurrentMSN_m(1,:) =  p.D1_mCurrentMSN_Q10^(.1*(37-23))*1e-4*(D1_mCurrentMSN_V_IC-p.D1_mCurrentMSN_vhalf)/(1-exp(-(D1_mCurrentMSN_V_IC-p.D1_mCurrentMSN_vhalf)/9))/(p.D1_mCurrentMSN_Q10^(.1*(37-23))*1e-4*(D1_mCurrentMSN_V_IC-p.D1_mCurrentMSN_vhalf)/(1-exp(-(D1_mCurrentMSN_V_IC-p.D1_mCurrentMSN_vhalf)/9))-p.D1_mCurrentMSN_Q10^(.1*(37-23))*1e-4*(D1_mCurrentMSN_V_IC-p.D1_mCurrentMSN_vhalf)/(1-exp((D1_mCurrentMSN_V_IC-p.D1_mCurrentMSN_vhalf)/9)))*randn(1,p.D1_Npop);
D2_V = zeros(nsamp,p.D2_Npop);
  D2_V(1,:) = zeros(1,p.D2_Npop);
D2_naCurrentMSN_m = zeros(nsamp,p.D2_Npop);
  D2_naCurrentMSN_m(1,:) =  0.32*(D2_naCurrentMSN_V_IC+54)/(1-exp(-(D2_naCurrentMSN_V_IC+54)/4))/(0.32*(D2_naCurrentMSN_V_IC+54)/(1-exp(-(D2_naCurrentMSN_V_IC+54)/4))+0.28*(D2_naCurrentMSN_V_IC+27)/(exp((D2_naCurrentMSN_V_IC+27)/5)-1))*randn(1,p.D2_Npop);
D2_naCurrentMSN_h = zeros(nsamp,p.D2_Npop);
  D2_naCurrentMSN_h(1,:) =  0.128*exp(-(D2_naCurrentMSN_V_IC+50)/18)/(0.128*exp(-(D2_naCurrentMSN_V_IC+50)/18)+4/(1+exp(-(D2_naCurrentMSN_V_IC+27)/5)))*randn(1,p.D2_Npop);
D2_kCurrentMSN_m = zeros(nsamp,p.D2_Npop);
  D2_kCurrentMSN_m(1,:) =  0.032*(D2_kCurrentMSN_V_IC+52)/(1-exp(-(D2_kCurrentMSN_V_IC+52)/5))/(0.032*(D2_kCurrentMSN_V_IC+52)/(1-exp(-(D2_kCurrentMSN_V_IC+52)/5))+0.5*exp(-(D2_kCurrentMSN_V_IC+57)/40))*randn(1,p.D2_Npop);
D2_mCurrentMSN_m = zeros(nsamp,p.D2_Npop);
  D2_mCurrentMSN_m(1,:) =  p.D2_mCurrentMSN_Q10^(.1*(37-23))*1e-4*(D2_mCurrentMSN_V_IC-p.D2_mCurrentMSN_vhalf)/(1-exp(-(D2_mCurrentMSN_V_IC-p.D2_mCurrentMSN_vhalf)/9))/(p.D2_mCurrentMSN_Q10^(.1*(37-23))*1e-4*(D2_mCurrentMSN_V_IC-p.D2_mCurrentMSN_vhalf)/(1-exp(-(D2_mCurrentMSN_V_IC-p.D2_mCurrentMSN_vhalf)/9))-p.D2_mCurrentMSN_Q10^(.1*(37-23))*1e-4*(D2_mCurrentMSN_V_IC-p.D2_mCurrentMSN_vhalf)/(1-exp((D2_mCurrentMSN_V_IC-p.D2_mCurrentMSN_vhalf)/9)))*randn(1,p.D2_Npop);
soma_soma_somaSomaiSYN_s = zeros(nsamp,p.soma_Npop);
  soma_soma_somaSomaiSYN_s(1,:) =  p.soma_soma_somaSomaiSYN_IC+p.soma_soma_somaSomaiSYN_IC_noise.*rand(1,p.soma_Npop);
D1_soma_somaMSNiSYN_s = zeros(nsamp,p.soma_Npop);
  D1_soma_somaMSNiSYN_s(1,:) =  p.D1_soma_somaMSNiSYN_IC+p.D1_soma_somaMSNiSYN_IC_noise.*rand(1,p.soma_Npop);
D2_soma_somaMSNiSYN_s = zeros(nsamp,p.soma_Npop);
  D2_soma_somaMSNiSYN_s(1,:) =  p.D2_soma_somaMSNiSYN_IC+p.D2_soma_somaMSNiSYN_IC_noise.*rand(1,p.soma_Npop);
D1_D1_gabaRecInputMSN_s = zeros(nsamp,p.D1_Npop);
  D1_D1_gabaRecInputMSN_s(1,:) = zeros(1,p.D1_Npop);
D2_D2_gabaRecInputMSN_s = zeros(nsamp,p.D2_Npop);
  D2_D2_gabaRecInputMSN_s(1,:) = zeros(1,p.D2_Npop);
% ###########################################################
% Numerical integration:
% ###########################################################
n=2;
for k=2:ntime
  t=T(k-1);
  soma_V_k1=p.soma_Iapp+((-(( p.soma_somaGolombK_gd.*(soma_somaGolombK_a(n-1,:).^3).*soma_somaGolombK_b(n-1,:).*(soma_V(n-1,:)-p.soma_somaGolombK_vk))))+((-(( p.soma_somaGolombKdr_gkdr.*(soma_somaGolombKdr_n(n-1,:).^p.soma_somaGolombKdr_nexp).*(soma_V(n-1,:)-p.soma_somaGolombKdr_vk))))+(((( p.soma_somaInput_tonic)))+((-(( p.soma_somaGolombNa_gna.*soma_somaGolombNa_h(n-1,:).*((( 1./(1+exp(-(soma_V(n-1,:)-p.soma_somaGolombNa_thetam)./11.5)))).^3).*(soma_V(n-1,:)-p.soma_somaGolombNa_vna))))+((-(( p.soma_somaLeak_gl*(soma_V(n-1,:)-p.soma_somaLeak_vl))))+((-(( (soma_soma_somaSomaiSYN_gsyn2.*(soma_soma_somaSomaiSYN_s(n-1,:)*soma_soma_somaSomaiSYN_mask).*(soma_V(n-1,:)-p.soma_soma_somaSomaiSYN_Esyn)))))+(((( p.soma_dend_iCOM_gCOM.*(dend_V(n-1,:)-soma_V(n-1,:))))))))))));
  soma_somaGolombK_a_k1= ((( 1./(1+exp(-(soma_V(n-1,:)+50)./20))))-soma_somaGolombK_a(n-1,:))./2;
  soma_somaGolombK_b_k1= ((( 1./(1+exp(-(soma_V(n-1,:)+70)./-6))))-soma_somaGolombK_b(n-1,:))./p.soma_somaGolombK_taub;
  soma_somaGolombKdr_n_k1= ((( 1./(1+exp(-(soma_V(n-1,:)+12.4)./6.8))))-soma_somaGolombKdr_n(n-1,:))./(( p.soma_somaGolombKdr_tau_mult*(0.087+11.4./(1+exp((soma_V(n-1,:)+14.6)./8.6))) .* (0.087+11.4./(1+exp(-(soma_V(n-1,:)-1.3)./18.7)))));
  soma_somaGolombNa_h_k1= ((( 1./(1+exp(-(soma_V(n-1,:)+58.3)./-6.7))))-soma_somaGolombNa_h(n-1,:))./(( 0.5 + 14./(1+exp(-(soma_V(n-1,:)+60)./-12))));
  dend_V_k1=p.dend_Iapp+((-(( dend_dendGolombK_gd_dend.*(dend_dendGolombK_a(n-1,:).^3).*dend_dendGolombK_b(n-1,:).*(dend_V(n-1,:)-p.dend_dendGolombK_vk))))+((-(( dend_dendGolombKdr_gkdr_dend.*(dend_dendGolombKdr_n(n-1,:).^p.dend_dendGolombKdr_nexp).*(dend_V(n-1,:)-p.dend_dendGolombKdr_vk))))+((-(( dend_dendGolombNa_gna_dend.*dend_dendGolombNa_h(n-1,:).*((( 1./(1+exp(-(dend_V(n-1,:)-p.dend_dendGolombNa_thetam)./11.5)))).^3).*(dend_V(n-1,:)-p.dend_dendGolombNa_vna))))+(((( dend_dendInput_tonic2)))+((-(( p.dend_dendLeak_gl*(dend_V(n-1,:)-p.dend_dendLeak_vl))))+((-(( (( (( p.dend_dendiMultiPoissonExp_g_esyn.*dend_dendiMultiPoissonExp_Ge(:, max(1,round(t/dt)))')).*(dend_V(n-1,:) - p.dend_dendiMultiPoissonExp_E_esyn))) + (( (( p.dend_dendiMultiPoissonExp_g_isyn.*dend_dendiMultiPoissonExp_Gi(:, max(1,round(t/dt)))')).*(dend_V(n-1,:) - p.dend_dendiMultiPoissonExp_E_isyn))))))+(((( p.dend_soma_iCOM_gCOM.*(soma_V(n-1,:)-dend_V(n-1,:)))))+(((( dend_dend_dendDendiGAP_g_GAP2.*sum((( (dend_V(n-1,:)'*ones(1,size(dend_V(n-1,:),2)))-(dend_V(n-1,:)'*ones(1,size(dend_V(n-1,:),2))))).*dend_dend_dendDendiGAP_mask,2)')))))))))));
  dend_dendGolombK_a_k1= ((( 1./(1+exp(-(dend_V(n-1,:)+50)./20))))-dend_dendGolombK_a(n-1,:))./2;
  dend_dendGolombK_b_k1= ((( 1./(1+exp(-(dend_V(n-1,:)+70)./-6))))-dend_dendGolombK_b(n-1,:))./p.dend_dendGolombK_taub;
  dend_dendGolombKdr_n_k1= ((( 1./(1+exp(-(dend_V(n-1,:)+12.4)./6.8))))-dend_dendGolombKdr_n(n-1,:))./(( p.dend_dendGolombKdr_tau_mult*(0.087+11.4./(1+exp((dend_V(n-1,:)+14.6)./8.6))) .* (0.087+11.4./(1+exp(-(dend_V(n-1,:)-1.3)./18.7)))));
  dend_dendGolombNa_h_k1= ((( 1./(1+exp(-(dend_V(n-1,:)+58.3)./-6.7))))-dend_dendGolombNa_h(n-1,:))./(( 0.5 + 14./(1+exp(-(dend_V(n-1,:)+60)./-12))));
  D1_V_k1=p.D1_Iapp+((-(( p.D1_naCurrentMSN_g_na*D1_naCurrentMSN_m(n-1,:).^3.*D1_naCurrentMSN_h(n-1,:).*(D1_V(n-1,:)-p.D1_naCurrentMSN_E_na))))+((-((p.D1_kCurrentMSN_g_k*D1_kCurrentMSN_m(n-1,:).^4.*(D1_V(n-1,:)-p.D1_kCurrentMSN_E_k))))+((-(( p.D1_mCurrentMSN_g_m*D1_mCurrentMSN_m(n-1,:).*(D1_V(n-1,:)-p.D1_mCurrentMSN_E_m))))+((-((p.D1_leakCurrentMSN_g_l*(D1_V(n-1,:)-p.D1_leakCurrentMSN_E_l))))+(((( p.D1_injectedCurrentD1_injectedCurrent + p.D1_injectedCurrentD1_sinmult*sin(2*pi*D1_injectedCurrentD1_freq*t/1000) + 0.1*p.D1_injectedCurrentD1_DA)))+(((( p.D1_noisyInputMSN_sigma_noise.*randn(1,p.D1_Npop).*sqrt(dt))))+((-(( (D1_soma_somaMSNiSYN_gsyn2.*(D1_soma_somaMSNiSYN_s(n-1,:)*D1_soma_somaMSNiSYN_mask).*(D1_V(n-1,:)-p.D1_soma_somaMSNiSYN_Esyn)))))+((-(( p.D1_D1_gabaRecInputMSN_g_gaba.*(D1_D1_gabaRecInputMSN_s(n-1,:)*D1_D1_gabaRecInputMSN_netcon).*(D1_V(n-1,:)-p.D1_D1_gabaRecInputMSN_E_gaba))))))))))));
  D1_naCurrentMSN_m_k1= (( 0.32*(D1_V(n-1,:)+54)./(1-exp(-(D1_V(n-1,:)+54)/4)))).*(1-D1_naCurrentMSN_m(n-1,:))-(( 0.28*(D1_V(n-1,:)+27)./(exp((D1_V(n-1,:)+27)/5)-1))).*D1_naCurrentMSN_m(n-1,:);
  D1_naCurrentMSN_h_k1= (( 0.128*exp(-(D1_V(n-1,:)+50)/18))).*(1-D1_naCurrentMSN_h(n-1,:))-(( 4./(1+exp(-(D1_V(n-1,:)+27)/5)))).*D1_naCurrentMSN_h(n-1,:);
  D1_kCurrentMSN_m_k1= (( 0.032*(D1_V(n-1,:)+52)./(1-exp(-(D1_V(n-1,:)+52)/5)))).*(1-D1_kCurrentMSN_m(n-1,:))-(( 0.5*exp(-(D1_V(n-1,:)+57)/40))).*D1_kCurrentMSN_m(n-1,:);
  D1_mCurrentMSN_m_k1= (( D1_mCurrentMSN_Qs*1e-4*(D1_V(n-1,:)-p.D1_mCurrentMSN_vhalf)./(1-exp(-(D1_V(n-1,:)-p.D1_mCurrentMSN_vhalf)/9)))).*(1-D1_mCurrentMSN_m(n-1,:)) - (( -D1_mCurrentMSN_Qs*1e-4*(D1_V(n-1,:)-p.D1_mCurrentMSN_vhalf)./(1-exp((D1_V(n-1,:)-p.D1_mCurrentMSN_vhalf)/9)))).*D1_mCurrentMSN_m(n-1,:);
  D2_V_k1=p.D2_Iapp+((-(( p.D2_naCurrentMSN_g_na*D2_naCurrentMSN_m(n-1,:).^3.*D2_naCurrentMSN_h(n-1,:).*(D2_V(n-1,:)-p.D2_naCurrentMSN_E_na))))+((-((p.D2_kCurrentMSN_g_k*D2_kCurrentMSN_m(n-1,:).^4.*(D2_V(n-1,:)-p.D2_kCurrentMSN_E_k))))+((-(( p.D2_mCurrentMSN_g_m*D2_mCurrentMSN_m(n-1,:).*(D2_V(n-1,:)-p.D2_mCurrentMSN_E_m))))+((-((p.D2_leakCurrentMSN_g_l*(D2_V(n-1,:)-p.D2_leakCurrentMSN_E_l))))+(((( p.D2_injectedCurrentD2_injectedCurrent + p.D2_injectedCurrentD2_sinmult*sin(2*pi*D2_injectedCurrentD2_freq*t/1000) - 0.1*p.D2_injectedCurrentD2_DA)))+(((( p.D2_noisyInputMSN_sigma_noise.*randn(1,p.D2_Npop).*sqrt(dt))))+((-(( (D2_soma_somaMSNiSYN_gsyn2.*(D2_soma_somaMSNiSYN_s(n-1,:)*D2_soma_somaMSNiSYN_mask).*(D2_V(n-1,:)-p.D2_soma_somaMSNiSYN_Esyn)))))+((-(( p.D2_D2_gabaRecInputMSN_g_gaba.*(D2_D2_gabaRecInputMSN_s(n-1,:)*D2_D2_gabaRecInputMSN_netcon).*(D2_V(n-1,:)-p.D2_D2_gabaRecInputMSN_E_gaba))))))))))));
  D2_naCurrentMSN_m_k1= (( 0.32*(D2_V(n-1,:)+54)./(1-exp(-(D2_V(n-1,:)+54)/4)))).*(1-D2_naCurrentMSN_m(n-1,:))-(( 0.28*(D2_V(n-1,:)+27)./(exp((D2_V(n-1,:)+27)/5)-1))).*D2_naCurrentMSN_m(n-1,:);
  D2_naCurrentMSN_h_k1= (( 0.128*exp(-(D2_V(n-1,:)+50)/18))).*(1-D2_naCurrentMSN_h(n-1,:))-(( 4./(1+exp(-(D2_V(n-1,:)+27)/5)))).*D2_naCurrentMSN_h(n-1,:);
  D2_kCurrentMSN_m_k1= (( 0.032*(D2_V(n-1,:)+52)./(1-exp(-(D2_V(n-1,:)+52)/5)))).*(1-D2_kCurrentMSN_m(n-1,:))-(( 0.5*exp(-(D2_V(n-1,:)+57)/40))).*D2_kCurrentMSN_m(n-1,:);
  D2_mCurrentMSN_m_k1= (( D2_mCurrentMSN_Qs*1e-4*(D2_V(n-1,:)-p.D2_mCurrentMSN_vhalf)./(1-exp(-(D2_V(n-1,:)-p.D2_mCurrentMSN_vhalf)/9)))).*(1-D2_mCurrentMSN_m(n-1,:)) - (( -D2_mCurrentMSN_Qs*1e-4*(D2_V(n-1,:)-p.D2_mCurrentMSN_vhalf)./(1-exp((D2_V(n-1,:)-p.D2_mCurrentMSN_vhalf)/9)))).*D2_mCurrentMSN_m(n-1,:);
  soma_soma_somaSomaiSYN_s_k1= -soma_soma_somaSomaiSYN_s(n-1,:)./p.soma_soma_somaSomaiSYN_tauD + ((1-soma_soma_somaSomaiSYN_s(n-1,:))/p.soma_soma_somaSomaiSYN_tauR).*(1+tanh(soma_V(n-1,:)/10));
  D1_soma_somaMSNiSYN_s_k1= -D1_soma_somaMSNiSYN_s(n-1,:)./p.D1_soma_somaMSNiSYN_tauD + ((1-D1_soma_somaMSNiSYN_s(n-1,:))/p.D1_soma_somaMSNiSYN_tauR).*(1+tanh(soma_V(n-1,:)/10));
  D2_soma_somaMSNiSYN_s_k1= -D2_soma_somaMSNiSYN_s(n-1,:)./p.D2_soma_somaMSNiSYN_tauD + ((1-D2_soma_somaMSNiSYN_s(n-1,:))/p.D2_soma_somaMSNiSYN_tauR).*(1+tanh(soma_V(n-1,:)/10));
  D1_D1_gabaRecInputMSN_s_k1= -D1_D1_gabaRecInputMSN_s(n-1,:)./p.D1_D1_gabaRecInputMSN_tau_gaba + 2*(1+tanh(D1_V(n-1,:)/4)).*(1-D1_D1_gabaRecInputMSN_s(n-1,:));
  D2_D2_gabaRecInputMSN_s_k1= -D2_D2_gabaRecInputMSN_s(n-1,:)./p.D2_D2_gabaRecInputMSN_tau_gaba + 2*(1+tanh(D2_V(n-1,:)/4)).*(1-D2_D2_gabaRecInputMSN_s(n-1,:));
  % ------------------------------------------------------------
  % Update state variables:
  % ------------------------------------------------------------
  soma_V(n,:)=soma_V(n-1,:)+dt*soma_V_k1;
  soma_somaGolombK_a(n,:)=soma_somaGolombK_a(n-1,:)+dt*soma_somaGolombK_a_k1;
  soma_somaGolombK_b(n,:)=soma_somaGolombK_b(n-1,:)+dt*soma_somaGolombK_b_k1;
  soma_somaGolombKdr_n(n,:)=soma_somaGolombKdr_n(n-1,:)+dt*soma_somaGolombKdr_n_k1;
  soma_somaGolombNa_h(n,:)=soma_somaGolombNa_h(n-1,:)+dt*soma_somaGolombNa_h_k1;
  dend_V(n,:)=dend_V(n-1,:)+dt*dend_V_k1;
  dend_dendGolombK_a(n,:)=dend_dendGolombK_a(n-1,:)+dt*dend_dendGolombK_a_k1;
  dend_dendGolombK_b(n,:)=dend_dendGolombK_b(n-1,:)+dt*dend_dendGolombK_b_k1;
  dend_dendGolombKdr_n(n,:)=dend_dendGolombKdr_n(n-1,:)+dt*dend_dendGolombKdr_n_k1;
  dend_dendGolombNa_h(n,:)=dend_dendGolombNa_h(n-1,:)+dt*dend_dendGolombNa_h_k1;
  D1_V(n,:)=D1_V(n-1,:)+dt*D1_V_k1;
  D1_naCurrentMSN_m(n,:)=D1_naCurrentMSN_m(n-1,:)+dt*D1_naCurrentMSN_m_k1;
  D1_naCurrentMSN_h(n,:)=D1_naCurrentMSN_h(n-1,:)+dt*D1_naCurrentMSN_h_k1;
  D1_kCurrentMSN_m(n,:)=D1_kCurrentMSN_m(n-1,:)+dt*D1_kCurrentMSN_m_k1;
  D1_mCurrentMSN_m(n,:)=D1_mCurrentMSN_m(n-1,:)+dt*D1_mCurrentMSN_m_k1;
  D2_V(n,:)=D2_V(n-1,:)+dt*D2_V_k1;
  D2_naCurrentMSN_m(n,:)=D2_naCurrentMSN_m(n-1,:)+dt*D2_naCurrentMSN_m_k1;
  D2_naCurrentMSN_h(n,:)=D2_naCurrentMSN_h(n-1,:)+dt*D2_naCurrentMSN_h_k1;
  D2_kCurrentMSN_m(n,:)=D2_kCurrentMSN_m(n-1,:)+dt*D2_kCurrentMSN_m_k1;
  D2_mCurrentMSN_m(n,:)=D2_mCurrentMSN_m(n-1,:)+dt*D2_mCurrentMSN_m_k1;
  soma_soma_somaSomaiSYN_s(n,:)=soma_soma_somaSomaiSYN_s(n-1,:)+dt*soma_soma_somaSomaiSYN_s_k1;
  D1_soma_somaMSNiSYN_s(n,:)=D1_soma_somaMSNiSYN_s(n-1,:)+dt*D1_soma_somaMSNiSYN_s_k1;
  D2_soma_somaMSNiSYN_s(n,:)=D2_soma_somaMSNiSYN_s(n-1,:)+dt*D2_soma_somaMSNiSYN_s_k1;
  D1_D1_gabaRecInputMSN_s(n,:)=D1_D1_gabaRecInputMSN_s(n-1,:)+dt*D1_D1_gabaRecInputMSN_s_k1;
  D2_D2_gabaRecInputMSN_s(n,:)=D2_D2_gabaRecInputMSN_s(n-1,:)+dt*D2_D2_gabaRecInputMSN_s_k1;
  n=n+1;
end
T=T(1:downsample_factor:ntime);

end


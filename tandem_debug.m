%% Test
T = 100;
h0 = 10.5;
delta_LN = 5;
TLs = 50;
x = rand(T,1);
sigma = zeros(T,1);

for i=2:T
   sigma(i) = sigma(i-1) + x(i)*(delta_LN*h0/TLs);
end

%plot(sigma);

%% Tandem
delta = 36;
T = 400;
h0 = 10.5;
tau = 0.25;
gamma = 0.5;
sigma_0 = -0.65;
delta_LN = 5;
deg_m = 0.0001;
omega = 10^5;%7.85*10^11;
TL_S = 141;
TL_G = 50;
PL_S = 40;
PL_G = 44;
NS = 50;
kf_max = 10^4;
kcat_max = 5.4*10^5;
k_w = 1;
k_l = 0.02;
k_r = 1;

sigma_tS = -0.65;
sigma_tG = -0.65;
sigma_pS = -0.65;
sigma_pG = -0.65;
mS = 0;
MG = 0;
EC_S = 0.1*10^-9;
EC_G = 0.1*10^-9;
k_seq = 1;
ECGECS = 0.1*10^-9;

R_tot = 10^-6 - EC_S + EC_G + ECGECS;
PLac_tot = 11*10^-9;% EC_S + ECGECS;
PTet_tot = 11*10^-9;% EC_G + ECGECS;

for i=1:T
    
% if(sigma_pS<sigma_0)Bt_pS = 1;
% else Bt_pS = 0;
% end
% if(sigma_pS>sigma_0)Bg_pS = 1;
% else Bg_pS = 0;
% end
% 
% if(sigma_tS<=sigma_0)Bt_tS = 1;
% else Bt_tS = 0;
% end
% if(sigma_tS>sigma_0)Bg_tS = 1;
% else Bg_tS = 0;
% end
% 
% if(sigma_pG<sigma_0)Bt_pG = 1;
% else Bt_pG = 0;
% end
% if(sigma_pG>sigma_0)Bg_pG = 1;
% else Bg_pG = 0;
% end
% 
% if(sigma_tG<sigma_0)Bt_tG = 1;
% else Bt_tG = 0;
% end
% if(sigma_tG>sigma_0)Bg_tG = 1;
% else Bg_tG = 0;
% end

if(sigma_tS(i)>sigma_0)BtS = -gamma;
else BtS = tau;
end

if(sigma_pG(i)>sigma_0)BpG = -gamma;
else BpG = tau;
end

if(sigma_pS(i)>sigma_0)BpS = -gamma;
else BpS = tau;
end

if(sigma_tG(i)>sigma_0)BtG = -gamma;
else BtG = tau;
end

%R_tot(i) = R + EC_S(i) + EC_G(i) + ECGECS(i);
PLac(i) = PLac_tot - EC_S(i) - ECGECS(i);
PTet(i) = PTet_tot - EC_G(i) - ECGECS(i);

sigma_tS(i+1) =  sigma_tS(i) + (delta)*(-(omega/2)*(kcat(sigma_tS(i), kcat_max, TL_S)*EC_S(i) ) + (h0/TL_S)*BtS );
sigma_pG(i+1) = sigma_pG(i) + delta*(-(omega/2)*(kf(sigma_pG(i), kf_max)*PTet(i)*(R_tot - (EC_S(i) + EC_G(i) + ECGECS(i)))) + (h0/PL_G)*BpG);

sigma_pS(i+1) = sigma_pS(i) + delta*((omega/2)*(kcat(sigma_pG(i), kcat_max, TL_G)*EC_G(i)*(TL_G/(2*(PL_S+NS))) - kf(sigma_pS(i), kf_max)*PLac(i)*(R_tot - (EC_S(i) + EC_G(i) + ECGECS(i)))-kcat(sigma_tS(i), kcat_max, TL_G)*EC_S(i)*(TL_G/(PL_S+NS))) + (h0/PL_S)*BpS);
sigma_tG(i+1) = sigma_tG(i) + delta*( -(omega/2)*(kcat(sigma_tS(i), kcat_max, TL_S)*EC_S(i)*(TL_S/(PL_S+NS+TL_G+TL_S)) + kcat(sigma_tG(i), kcat_max, TL_G)*EC_G(i) + kf(sigma_pS(i), kf_max)*PLac(i)*(R_tot - (EC_S(i) + EC_G(i) + ECGECS(i)))*(PL_S/(2*(TL_G+NS)))) + (h0/TL_S)*BtG);

mS(i+1) = mS(i) + (delta)*(kcat(sigma_tS(i), kcat_max, TL_S)*EC_S(i)+k_w*ECGECS(i)-deg_m*mS(i));
MG(i+1) = MG(i) + delta*(kcat(sigma_tG(i), kcat_max, TL_S)*EC_G(i)+k_w*ECGECS(i)-deg_m*MG(i));
EC_S(i+1) = EC_S(i) ;%+ (delta/1)*(kf(sigma_pS(i), kf_max)*(R_tot - (EC_S(i) + EC_G(i) + ECGECS(i)))*PLac(i)-(k_r+kcat(sigma_tS(i), kcat_max, TL_S))*EC_S(i));
EC_G(i+1) = EC_G(i) ;%+ (delta/1)*(kf(sigma_pG(i), kf_max)*(R_tot - (EC_S(i) + EC_G(i) + ECGECS(i)))*PTet(i)-(k_r+kcat(sigma_tG(i), kcat_max, TL_G)+kseq(sigma_tG(i))+k_l)*EC_G(i));
ECGECS(i+1) = ECGECS(i) ;%+ (delta/1)*(k_l*EC_G(i) - k_w*ECGECS(i));

X(i) = EC_S(i) + EC_G(i) + ECGECS(i);
Y(i) = kcat(sigma_tS(i), kcat_max, TL_S);

end
figure;
hold on;
%plot(Y);
%plot(mS);
plot(MG);
% plot(sigma_tS);
% plot(sigma_pG);
% plot(sigma_pS);
% plot(sigma_tG);

%plot(EC_S);
%plot(EC_G);
%plot(ECGECS);

legend('show')

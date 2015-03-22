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
delta = 0.001;
T = 10000;
h0 = 10.5;
tau = 0.25;
gamma = 0.5;
sigma_0 = -0.65;
delta_LN = 5;
delta_m = 0.05;
omega = 100000/1;
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

R = 0.1;
R_tot = 0.1;
PLac = 0.1;
PTet = 0.1;

sigma_tS = -0.65;
sigma_tG = -0.65;
sigma_pS = -0.65;
sigma_pG = -0.65;
mS = 0;
MG = 0.1;
EC_S = 0.1;
EC_G = 0.1;
k_seq = 0.1;
ECGECS = 0.1;

for i=1:T
    
if(sigma_pS<sigma_0)Bt_pS = 1;
else Bt_pS = 0;
end
if(sigma_pS>sigma_0)Bg_pS = 1;
else Bg_pS = 0;
end

if(sigma_tS<sigma_0)Bt_tS = 1;
else Bt_tS = 0;
end
if(sigma_tS>sigma_0)Bg_tS = 1;
else Bg_tS = 0;
end

if(sigma_pG<sigma_0)Bt_pG = 1;
else Bt_pG = 0;
end
if(sigma_pG>sigma_0)Bg_pG = 1;
else Bg_pG = 0;
end

if(sigma_tG<sigma_0)Bt_tG = 1;
else Bt_tG = 0;
end
if(sigma_tG>sigma_0)Bg_tG = 1;
else Bg_tG = 0;
end

sigma_tS(i+1) = sigma_tS(i) + delta*(-(omega/2)*(kcat(sigma_tS(i), kcat_max, TL_S)*EC_G(i) + (h0/TL_S)*(tau*Bt_tS-gamma*Bg_tS)));
sigma_pG(i+1) = sigma_pG(i) + delta*(-(omega/2)*(kf(sigma_pG(i), kf_max)*PTet*R + (h0/PL_G)*(tau*Bt_pG-gamma*Bg_pG)));

sigma_pS(i+1) = sigma_pS(i) + delta*((omega/2)*(kcat(sigma_pG(i), kcat_max, TL_G)*EC_G(i)*(TL_G/(2*(PL_S+NS))) - kf(sigma_pS(i), kf_max)*PLac*R-kcat(sigma_tS(i), kcat_max, TL_G)*EC_S(i)*(TL_G/(PL_S+NS))) + (h0/PL_S)*(tau*Bt_pS-gamma*Bg_pS));
sigma_tG(i+1) = sigma_tG(i) + delta*(-(omega/2)*(kcat(sigma_tG(i), kcat_max, TL_G)*EC_G(i) + kcat(sigma_tG(i), kcat_max, TL_G)*EC_S(i)*(TL_S/(PL_S+NS+TL_G+TL_S)) + kf(sigma_pS(i), kf_max)*PLac*R*(PL_S/(2*(TL_G+NS))) + (h0/TL_S)*(tau*Bt_tG-gamma*Bg_tG)));

mS(i+1) = mS(i) + delta*(kcat(sigma_tS(i), kcat_max, TL_S)*EC_S(i)+k_w*ECGECS(i)-delta_m*mS(i));
MG(i+1) = MG(i) + delta*(kcat(sigma_tG(i), kcat_max, TL_S)*EC_G(i)+k_w*ECGECS(i)-delta_m*mS(i));
EC_S(i+1) = EC_S(i) + delta*(kf(sigma_pS(i), kf_max)*R*(PLac-EC_S(i)-ECGECS(i))-(k_r+kcat(sigma_tS(i), kcat_max, TL_S))*EC_S(i));
EC_G(i+1) = EC_G(i) + delta*(kf(sigma_pG(i), kf_max)*R*(PTet-EC_G(i)-ECGECS(i))-(k_r+kcat(sigma_tG(i), kcat_max, TL_G)+k_seq+k_l)*EC_G(i));
ECGECS(i+1) = ECGECS(i) + delta*(k_l*EC_G(i) - k_w*ECGECS(i));

R_tot = R + EC_S(i) + EC_G(i) + ECGECS(i);
end
figure;
hold on;
% plot(sigma_pS);
% plot(sigma_tS);
% plot(sigma_pG);
% plot(sigma_tG);
plot(MG);
plot(mS);


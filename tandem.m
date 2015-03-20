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
h0 = 10.5;
tau = 0.25;
gamma = 0.5;
sigma_0 = -0.65;
delta_LN = 5;
delta_m = 0.05;
omega = 100/1;
TL_S = 141;
TL_G = 50;
PL_S = 40;
PL_G = 44;
NS = 50;
kf_max = 10^4;
kcat_max = 5.4*10^5;

R = 1;
PLac = 1;
PTet = 1;

sigma_tS = 1;
sigma_tG = 1;
sigma_pS = 1;
sigma_pG = 1;
mS = 1;
MG = 1;
EC_S = 1;
EC_G = 1;
GC = 1;

for i=1:10
if(sigma_pS<sigma_0)Bt = 1;
else Bt = 0;
end
if(sigma_pS>sigma_0)Bg = 1;
else Bg = 0;
end
sigma_pS = (omega/2)*(kcat(sigma_pG, kcat_max, TL_G)*EC_G*(TL_G/(2*(PL_S+NS)))-kf(sigma_pS, kf_max)*PLac*R-kcat(sigma_tS, kcat_max, TL_G)*EC_S*(TL_G/(PL_S+NS)))+h0/PL_S*(tau*Bt-gamma*Bg)
end



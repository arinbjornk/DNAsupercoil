% Supercoiling model
% Arinbjorn Kolbeinsson
% Imperial College London
% 2015

% Equations adapted from "Modeling the Effects of Compositional Context on Promoter Activity in an E. coli Extract based Transcription-Translation System" Enoch Yeung, Andrew Ng, Jongmin Kim,
% Zachary Z. Sun, and Richard M. Murray.
% http://www.cds.caltech.edu/~murray/papers/yeu+14-cdc.html

%% Parameters
delta = 0.036; %step size (seconds)
T = 400000; %number of cycles
h0 = 10.5;
tau = 0.25;
gamma = 0.5;
sigma_0 = -0.65;
deg_m = 0.0001;
omega = 7.85*10^11;
TL_S = 681;
TL_G = 720; 
PL_S = 101;
PL_G = 1210;
NS = 105;
kf_max = 0.2*10^-5;
kcat_max = 5.4*10^-4;
k_w = 10^-9;
k_l = 0.02;
k_r = 0.01;

%Initial conditions
sigma_tS = -0.65;
sigma_tG = -0.65;
sigma_pS = -0.65;
sigma_pG = -0.65;
mS = 0;
MG = 0;
k_seq = 1;
EC_S = 0;
EC_G = 0;
ECGECS = 0;

R_tot = 10^-6 - EC_S + EC_G + ECGECS;
PLac_tot = 11*10^-9;
PTet_tot = 11*10^-9;

%Rate dynamics
for i=1:T

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


PLac(i) = PLac_tot - EC_S(i) - ECGECS(i);
PTet(i) = PTet_tot - EC_G(i) - ECGECS(i);

sigma_tS(i+1) =  sigma_tS(i) + (delta)*(-(omega/2)*(kcat(sigma_tS(i), kcat_max, TL_S)*EC_S(i) ) + (h0/TL_S)*BtS );
sigma_pG(i+1) = sigma_pG(i) + delta*(-(omega/2)*(kf(sigma_pG(i), kf_max)*PTet(i)*(R_tot - (EC_S(i) + EC_G(i) + ECGECS(i)))) + (h0/PL_G)*BpG);

sigma_pS(i+1) = sigma_pS(i) + delta*((omega/2)*(kcat(sigma_pG(i), kcat_max, TL_G)*EC_G(i)*(TL_G/(2*(PL_S+NS))) - kf(sigma_pS(i), kf_max)*PLac(i)*(R_tot - (EC_S(i) + EC_G(i) + ECGECS(i)))-kcat(sigma_tS(i), kcat_max, TL_G)*EC_S(i)*(TL_G/(PL_S+NS))) + (h0/PL_S)*BpS);
sigma_tG(i+1) = sigma_tG(i) + delta*( -(omega/2)*(kcat(sigma_tS(i), kcat_max, TL_S)*EC_S(i)*(TL_S/(PL_S+NS+TL_G+TL_S)) + kcat(sigma_tG(i), kcat_max, TL_G)*EC_G(i) + kf(sigma_pS(i), kf_max)*PLac(i)*(R_tot - (EC_S(i) + EC_G(i) + ECGECS(i)))*(PL_S/(2*(TL_G+NS)))) + (h0/TL_S)*BtG);

mS(i+1) = mS(i) + (delta)*(kcat(sigma_tS(i), kcat_max, TL_S)*EC_S(i)+k_w*ECGECS(i)-deg_m*mS(i));
MG(i+1) = MG(i) + delta*(kcat(sigma_tG(i), kcat_max, TL_S)*EC_G(i)+k_w*ECGECS(i)-deg_m*MG(i));
EC_S(i+1) = EC_S(i) + (delta)*(kf(sigma_pS(i), kf_max)*(R_tot - (EC_S(i) + EC_G(i) + ECGECS(i)))*PLac(i)-(k_r+kcat(sigma_tS(i), kcat_max, TL_S))*EC_S(i));
EC_G(i+1) = EC_G(i) + (delta)*(kf(sigma_pG(i), kf_max)*(R_tot - (EC_S(i) + EC_G(i) + ECGECS(i)))*PTet(i)-(k_r+kcat(sigma_tG(i), kcat_max, TL_G)+kseq(sigma_tG(i))+k_l)*EC_G(i));
ECGECS(i+1) = ECGECS(i) + (delta)*(k_l*EC_G(i) - k_w*ECGECS(i));

end

figure;
%subplot(2,2,1);
hold on;
plot(mS,'LineWidth',4,'Color',[0 0.85 0.5]);
plot(MG,'LineWidth',4,'Color','r');
    xlab = 0:0.5:4.5;
    ylim([0 10^-19])
    xp = xlab*12;
%     set(gca,'XTick',xp); % Change x-axis ticks
     set(gca,'XTickLabel',xlab); % Change x-axis ticks labels to desired values.
    xlabel('Time/h') % x-axis label
    ylabel('Protein concentration /nM') % y-axis label
    set(gca,'FontSize',18,'FontWeight','bold')



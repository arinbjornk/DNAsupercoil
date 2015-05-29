TL_U = 141;
TL_D = 50;
s = [-800:-TL_U];
g = [-TL_U:0];
s2 = [0:800];
T = zeros(length(g));

y = TL_U./(2*(TL_D+s))-0.65;
y2 = TL_U./(2*(TL_D+s2))-0.65;
plot(s,y);
hold on;
plot(s2,y2);
plot(g,T,'LineWidth',3);
TL_U = 201;
TL_D = 50;
s = [-800:-TL_U];
g = [-TL_U:0];
s2 = [0:800];
s3 = [200:800];
T = ones(length(g))*-0.65;
dna = ones(1601)*-0.65;
dnax = [-800:800];
ins = ones(31)*-0.65;
insx = [185:215];

% insulator
for i=1:600
u(i+1) = sin(10*i);
    if(u(i+1) > 0) alpha = 1;
    else alpha = 0;
    end
end

%plotting
y = TL_U./(2*(TL_D+s))-0.65;
y2 = TL_U./(2*(TL_D+s2))-0.65;
y3 = (y2(200:800)+0.65)/4-0.65;
figure;
hold on;
plot(dnax,dna,'color','k','LineWidth',2);
plot(g,T,'LineWidth',10,'color',[0.7,0.4,0]);
plot(insx,ins,'LineWidth',10,'color','b');
plot(s,y,'LineWidth',3,'color',[0,0.7,0]);
plot(s2(1:201),y2(1:201),'LineWidth',3,'color',[0,0.7,0]);
plot(s3,y3,'LineWidth',3,'color',[0,0.6,1]);
    xlabel('Distance /bp') % x-axis label
    ylabel('Supercoiling density \sigma') % y-axis label
    set(gca,'FontSize',18,'FontWeight','bold')

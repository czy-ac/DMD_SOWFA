function [meandeltainstp]=evaluatetimevaryingerror(X,statesrebuild,dirdmd)

timevecc=size(X,2);
timevec=1:1:timevecc;

%conversion of time instant to continuous time in simulation
for l=1:length(timevec)
    k=timevec(l)+250;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    timeveccont(l)=integ;
end


for t=1:length(timevec)
    deltastate(:,t)=abs(X(:,t)-real(statesrebuild(:,t)));
    deltastatep(:,t)=abs(X(:,t)-real(statesrebuild(:,t)))./abs(X(:,t));
    meandeltainst(:,t)=mean(deltastate(:,t));
    meandeltainstp(:,t)=mean(deltastatep(:,t))*100;
    %deltainst(:,t)=norm(deltastate(:,t));
end

figure(550)

sid=scatter(timevec,meandeltainstp,'o');
hold on
sid.MarkerFaceColor = [0.2 0.6 0.8];
sid.MarkerEdgeColor = [0.2 0.6 0.8];
pid=plot(meandeltainstp,'LineWidth',1.1','color','blue');
pid.LineStyle='- -';
pid.Color=[0.8 0.8 1];

xlabel('Time [minutes]');
axis([-0.3667 750 0 max(meandeltainstp)+5])
ax=gca;
set(gca,'XTick',(0-0.3667:60:max(timevec))+1)
ax.XTickLabel = {'8','10','12','14','16','18','20','22','24','26','28','30','32','34'};
ylabel('Mean instantaneous deviaton from true state (%) ');
title('Dynamic Mode Decomposition state reconstruction: mean instataneous error with relation to true state');
set(gca, 'FontSize', 14);
grid on
grid minor
hold off
ytickformat('percentage')

set(gcf,'color','w','Position', get(0, 'Screensize'));    
export_fig(figure(550),strcat(dirdmd,'/image','statestimeerror'),'-nocrop','-m2'); 

function []=evauatemodelerror(states, statesrebuild,D,dirdmd,filename,dirName,x,y,z,Decimate)

 n=1;
    cases = dirName;

    [nTurbine,time4,dt,nVal,azi]           = readTurbineOutputGlobal(cases{n},'rotorAzimuth');
    [nTurbine,time5,dt,nVal,rotorPower{n}] = readTurbineOutputGlobal(cases{n},'rotorPower');
    [nTurbine,time6,dt,nVal,yawangle]      = readTurbineOutputGlobal(cases{n},'nacelleYaw');
    
    azirs(:,:)          =resample(azi(end-750*10:1:end,1:end),1,10);
    rotorPowerrs{n}(:,:)=resample(rotorPower{n}(end-750*10:1:end,1:end),1,10);
    yawanglers(:,:)     =resample(yawangle(end-750*10:1:end,1:end),1,10);

    load(filename);
    Uups=9; %[m/s]
    
    [xx,yy,zz]=resamplegrid(x,y,z, Decimate);
    [Xm_sh,Ym_sh,Zm_sh] = meshgrid(xx-500,(yy-500),zz);
    X = length(xx);
    Y = length(yy);
    Z = length(zz);
    
    %% 1. EVALUATE TIME AVERAGED ERROR
    delta=real(states-real(statesrebuild));
    
    for i=1:size(delta,1)
        taerror(i)=mean(delta(i,:));
    end
    taerror=taerror';
    
    fig503= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
    UmeanAbs_sh_u = reshape(taerror,Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Usq=squeeze(Usecu);
    plotturbinefromabove(0, 0, 0, D);
    hold on
    plotturbinefromabove(0, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Usq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Usq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([4 1 1])
    titlee=title(['Time-averaged error']);
    titlee.FontSize=18;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = ' (u_{SOWFA} - u_{DMD} ) / U ';
    set(gca,'fontsize', 16) 
    hold off 
    
    set(gcf,'color','w','Position', get(0, 'Screensize'));    
    shg
    export_fig(fig503,strcat(dirdmd,'/image','onetimeerror'),'-nocrop','-m2'); 
    
    %% EVALUATE INSTANTANEOUS ERROR
    %delta=abs((real(states(:,1:end-1))-real(statesrebuild)))./abs((states(:,1:end-1)))*100;
    
    %% First figure
    fig504= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
    
    i=405;
    subplot(5,2,1)
    plotsnapshothhdeviation(delta,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=407;
    subplot(5,2,2)
    plotsnapshothhdeviation(delta,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=420;
    subplot(5,2,3)
    plotsnapshothhdeviation(delta,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=430;
    subplot(5,2,4)
    plotsnapshothhdeviation(delta,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=440;
    subplot(5,2,5)
    plotsnapshothhdeviation(delta,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=450;
    subplot(5,2,6)
    plotsnapshothhdeviation(delta,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=460;
    subplot(5,2,7)
    plotsnapshothhdeviation(delta,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=470;
    subplot(5,2,8)
    plotsnapshothhdeviation(delta,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=480;
    subplot(5,2,9)
    plotsnapshothhdeviation(delta,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=490;
    subplot(5,2,10)
    plotsnapshothhdeviation(delta,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    titlee=suptitle(['Instantaneous deviation of DMD flow field reconstruction with relation to SOWFA for different time instances']);
    titlee.FontSize=18;
    
    set(gcf,'color','w','Position', get(0, 'Screensize'));    
    shg
    export_fig(fig504,strcat(dirdmd,'/image','errorvariationtime'),'-nocrop','-m2'); 
    
    %% INSTANTANEOUS SHOT DURING CRITICAL MOMENT
    fig505= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
    
    i=470;
    subplot(3,1,1)
    plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    subplot(3,1,2)
    plotsnapshothh(statesrebuild,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    subplot(3,1,3)
    plotsnapshothhdeviation(delta,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
     k=i+250;
     number=k/30;
     integ=floor(number);
     fract=number-integ;
     
     minutos=integ;
     segundos=60*fract;
     
    titlee=suptitle(['Instantaneous deviation of DMD flow field reconstruction with relation to SOWFA. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=18;
    
    set(gcf,'color','w','Position', get(0, 'Screensize'));    
    shg
    export_fig(fig505,strcat(dirdmd,'/image','errorvariationtime2'),'-nocrop','-m2'); 
    
    
    
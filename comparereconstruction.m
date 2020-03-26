function []=comparereconstruction(states, statesrebuild,D,dirdmd,x,y,z,Decimate,dirName)

    n=1;
    cases = dirName;

    [nTurbine,time4,dt,nVal,azi]           = readTurbineOutputGlobal(cases{n},'rotorAzimuth');
    [nTurbine,time5,dt,nVal,rotorPower{n}] = readTurbineOutputGlobal(cases{n},'rotorPower');
    [nTurbine,time6,dt,nVal,yawangle]      = readTurbineOutputGlobal(cases{n},'nacelleYaw');
    
    azirs(:,:)          =resample(azi(end-750*10:1:end,1:end),1,10);
    rotorPowerrs{n}(:,:)=resample(rotorPower{n}(end-750*10:1:end,1:end),1,10);
    yawanglers(:,:)     =resample(yawangle(end-750*10:1:end,1:end),1,10);

    Uups=9; %[m/s]
    
    [xx,yy,zz]=resamplegrid(x,y,z, Decimate);
    [Xm_sh,Ym_sh,Zm_sh] = meshgrid(xx-500,(yy-500),zz);
    X = length(xx);
    Y = length(yy);
    Z = length(zz);
    
    %% First figure
    fig500= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    
    i=405;
    subplot(5,2,1)
    plotsnapshothh(statesrebuild,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    subplot(5,2,2)
    plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=407;
    subplot(5,2,3)
    plotsnapshothh(statesrebuild,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    subplot(5,2,4)
    plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=420;
    subplot(5,2,5)
    plotsnapshothh(statesrebuild,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    subplot(5,2,6)
    plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=430;
    subplot(5,2,7)
    plotsnapshothh(statesrebuild,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    subplot(5,2,8)
    plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=440;
    subplot(5,2,9)
    plotsnapshothh(statesrebuild,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    subplot(5,2,10)
    plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    titlee=suptitle(['DMD flow field reconstruction                                                                                SOWFA flow field']);
    titlee.FontSize=18;
    
    set(gcf,'color','w','Position', get(0, 'Screensize'));    
    shg
    export_fig(fig500,strcat(dirdmd,'/image','reconsversusSOWFA'),'-nocrop','-m2'); 
 
   
    %% SECOND FIGURE
    fig502= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    
    i=450;
    subplot(5,2,1)
    plotsnapshothh(statesrebuild,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    subplot(5,2,2)
    plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=460;
    subplot(5,2,3)
    plotsnapshothh(statesrebuild,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    subplot(5,2,4)
    plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=470;
    subplot(5,2,5)
    plotsnapshothh(statesrebuild,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    subplot(5,2,6)
    plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=480;
    subplot(5,2,7)
    plotsnapshothh(statesrebuild,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    subplot(5,2,8)
    plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    i=490;
    subplot(5,2,9)
    plotsnapshothh(statesrebuild,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    subplot(5,2,10)
    plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)
    
    
     
    titlee=suptitle(['DMD flow field reconstruction                                                                                 SOWFA flow field']);
    titlee.FontSize=18;
    
    
    set(gcf,'color','w','Position', get(0, 'Screensize'));    
    shg
    export_fig(fig502,strcat(dirdmd,'/image','reconsversusSOWFA2'),'-nocrop','-m2'); 
    
    
    
function [FITje,OMEGA,DAMPING,fig1]=evaluatemodel(sys_red,si,Inputs, Outputs,FITje,OMEGA,DAMPING,purpose)
dN=2;
%Estimate the initial state, given the estimated system matrices, and a  set of input/output data.

a=strcmp(purpose, 'identification');

if a
    [xo]=dinit(sys_red{si}.A,sys_red{si}.B,sys_red{si}.C,sys_red{si}.D,[Inputs]',[Outputs]');

    %% VARIANCE ACCOUNTED FOR IN THE MODEL FOR IDENTIFICATION DATA
    %time response of dynamic system based on (1) derived state space
    %model (2) system inputs (3) computed initial conditions 
    ysim=lsim(sys_red{si}, [Inputs]',[],xo);   
    maxL(si)=max(abs(eig(sys_red{si}.A)));

    %compute the variance accounted for in the model based on simualton
    %of model and real output (for each time instant
    FITje(:,si)=vaf(ysim,[Outputs]');

    %% GRAPHICAL VISUALISAITON OF MODEL PREDICTION AND TRUE SIMULATION (IDENTIFICATION) RESULTS
    fig1=figure(1000+si);
    fig1.Visible='off';
    set(gcf,'color','w','Position', get(0, 'Screensize'));
    subplot(2,1,1)
    plot([Outputs(1,1:end-1);]','LineWidth',1.6'); %1
    hold on; 
    plot(ysim(:,1),'g--','LineWidth',1.6) %3
    grid on
    xlabel('Time instant [ ]')
    ylabel(' \Omega_1 [rad/s]')
    title(['Model fitness: rotor speed for first turbine. VAF of ',num2str(FITje(1,si)),' % '])
    legend({'Real simulated rotor speed','Model Output'},'Location','southeast')
    set(gca,'fontsize', 14)

    subplot(2,1,2)
    plot([Outputs(2,1:end-1)]','LineWidth',1.6'); %1
    hold on; 
    plot(ysim(:,2),'g--','LineWidth',1.6) %2
    grid on
    xlabel('Time instant [ ]')
    ylabel(' \Omega_2 [rad/s]')
    title(['Model fitness: rotor speed for second turbine. VAF of ',num2str(FITje(2,si)),' % '])
    legend({'Real simulated rotor speed','Model Output'},'Location','southeast') 
    set(gca,'fontsize', 14)

    %% GRAPHICAL VISUALISAITON OF MODEL PREDICTION AND TRUE SIMULATION (VALIDATION) RESULTS
else
    
     [xo_val]=dinit(sys_red{si}.A,sys_red{si}.B,sys_red{si}.C,sys_red{si}.D,[Inputs]',[Outputs]');
     ysim_val=lsim(sys_red{si}, [Inputs]',[],xo_val);  
     FITje(:,si)=vaf(ysim_val,[Outputs]');  
    
    fig1=figure(2000+si);
    fig1.Visible='off';
    set(gcf,'color','w','Position', get(0, 'Screensize'));
    subplot(2,1,1)
    plot([Outputs(1,1:end-1);]','LineWidth',1.6'); %1
    hold on; 
    plot(ysim_val(:,1),'g--','LineWidth',1.6) %3
    grid on
    xlabel('Time instant [ ]')
    ylabel(' \Omega_1 [rad/s]')
    title(['Model fitness: rotor speed for first turbine. VAF of ',num2str(FITje(1,si)),' % '])
    legend({'Real simulated rotor speed','Model Output',''},'Location','southeast')
    set(gca,'fontsize', 14)
             
    subplot(2,1,2,'Visible','off')
    plot([Outputs(2,1:end-1)]','LineWidth',1.6'); %1
    hold on; 
    plot(ysim_val(:,2),'g--','LineWidth',1.6)
    grid on
    xlabel('Time instant [ ]')
    ylabel(' \Omega_2 [rad/s]')
    title(['Model fitness: rotor speed for second turbine. VAF of ',num2str(FITje(2,si)),' % '])
    legend({'Real simulated rotor speed','Model Output',''},'Location','southeast') 
    set(gca,'fontsize', 14)

end


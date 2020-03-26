function [yf1, yf2]=evaluatepredictionpower(sys_red, modeltouse,Inputs_val, Outputs_val,dT,pasttime,stepsahead,forecasthorizon)

model=sys_red{modeltouse};
timevec=1:1:length(Inputs_val);
idfrw=iddata( Outputs_val',Inputs_val', dT);

%% PREDICTION: to evaluate if our model is a good predictor model
%yp = predict(sys,data,K) predicts the output of an identified model sys,
...K steps ahead using the measured input-output data.
time1=1;
time2=751;
inputdata=Inputs_val(time1:time2);
outputdata=Outputs_val(:,time1:time2);
data = iddata(outputdata',inputdata',dT);
warning off
[ypred xopred sys_pred]=predict(model,data,stepsahead);
warning on
y1=ypred.OutputData(:,1);
y2=ypred.OutputData(:,2);

%% LINEAR INSTANT SIMULATION
[xo]=dinit(model.A,model.B,model.C,model.D,Inputs_val',Outputs_val');
ysim=lsim(model, Inputs_val',[],xo);   

%% FORECASTING
%yf = forecast(sys,PastData,K) forecasts the output of an identified time series model 
...sys, K steps into the future using past measured data, PastData.
...forecast performs prediction into the future, in a time range beyond 
...the last instant of measured data. 
...In contrast, the predict command predicts the response of an identified model...
...over the time span of measured data. Use predict to determine if the predicted 
...result matches the observed response of an estimated model. 
...If sys is a good prediction model, consider using it with forecast.
time3=pasttime(1);
time4=pasttime(2);
%inputdatapast=idfrw.InputData(time3:time4);
%outputdatapast=idfrw.OutputData(time3:time4,:);
pastdata = idfrw(time3:time4);
%pastdata = iddata(outputdatapast,inputdatapast',dT);
sys=idss(model);
futureinput=[idfrw.InputData(time4+1:time4+forecasthorizon)];
%futureinputs=iddata(futureinput,[],2);
[xo]=dinit(sys_red{modeltouse}.A,sys_red{modeltouse}.B,sys_red{modeltouse}.C,...
    sys_red{modeltouse}.D,[Inputs_val(time4+1:end)]',[Outputs_val(:,time4+1:end)]');
opt = forecastOptions('InitialCOndition',idpar(xo));
yfuture = forecast(sys,pastdata,forecasthorizon,futureinput);
yf1=yfuture.OutputData(:,1);
yf2=yfuture.OutputData(:,2);
uf1=yfuture.InputData(:,1);

%% Plot de initial data 
figure(600)
set(gcf,'color','w','Position', get(0, 'Screensize'));

% subplot(3,1,1)
% plot(Inputs_val(1,1:end-1)','LineWidth',1.6'); %1
% hold on; 
% %plot(ysim(:,1),'g--','LineWidth',1.6) %3
% grid on
% ylabel(' \gamma_1 [deg]')
% %title(['Model fitness: rotor speed for first turbine. VAF of ',num2str(FITje(1,si)),' % '])
% %legend({'Real simulated rotor speed','Model Output'},'Location','southeast')
% set(gca,'fontsize', 14)
% l=line([time4 time4] , [min(Inputs_val(1,1:end-1))*(1.2) max(Inputs_val(1,1:end-1))*(1.2)]);
% l.Color=[0 0 0];
% l.LineWidth=2.5;
% l.LineStyle='- -';
% axis([0 length(Inputs_val(1:1:end-1)) min(Inputs_val(1,1:end-1))*(1.2) max(Inputs_val(1,1:end-1))*(1.2)]);
% f=plot(time4+1:1:time4+forecasthorizon, uf1);
% f.Color='y';
% f.Color=[1 0.5 0 ];
% f.LineWidth=1.6;
% f.LineStyle='-';

subplot(2,1,1)
plot(Outputs_val(1,1:end-1)','LineWidth',1.6'); %1
hold on; 
p=plot(ysim(:,1),'g--');
p.Color=[0.2 0.7 0.2];
P.LineWidth=1.6;
grid on
ylabel(' \Omega_1 [rad/s]')
title(['\Omega_1 forecasting from instant ',num2str(time4),' and ',num2str(forecasthorizon),' steps ahead in the future'])
set(gca,'fontsize', 14)
plot(y1,'r--');
l=line([time4 time4] , [min(Outputs_val(1,1:end-1))*(1.2) max(Outputs_val(1,1:end-1))*(1.2)]);
l.Color=[0 0 0];
l.LineWidth=2.5;
l.LineStyle='- -';
axis([0 length(Inputs_val(1:1:end-1)) min(Outputs_val(1,1:end-1))*(1.2) max(Outputs_val(1,1:end-1))*(1.2)]);
f=plot(time4+1:1:time4+forecasthorizon, yf1);
f.Color=[1 0.5 0 ];
f.LineWidth=1.6;
f.LineStyle='-';
legend({'\Omega_1 SOWFA','\Omega_1 simulated','\Omega_1 predicted','forecast frontier','\Omega_1 forecasted'},'Location','bestoutside','Orientation','horizontal')
legend('boxoff')

subplot(2,1,2)
plot(Outputs_val(2,1:end-1)','LineWidth',1.6'); %1
hold on; 
p=plot(ysim(:,2),'g--') ;
p.Color=[0.2 0.7 0.2];
p.LineWidth=1.6;
grid on
ylabel(' \Omega_2 [rad/s]')
title(['\Omega_2 forecasting from instant ',num2str(time4),' and ',num2str(forecasthorizon),' steps ahead in the future'])
set(gca,'fontsize', 14)
plot(y2,'r--');
l=line([time4 time4] , [min(Outputs_val(2,1:end-1))*(1.2) max(Outputs_val(2,1:end-1))*(1.2)]);
l.Color=[0 0 0];
l.LineWidth=2.5;
l.LineStyle='- -';
axis([0 length(Inputs_val(1:1:end-1)) min(Outputs_val(2,1:end-1))*(1.2) max(Outputs_val(2,1:end-1))*(1.2)]);
f=plot(timevec(time4+1:1:time4+forecasthorizon), yf2);
f.Color='y';
f.Color=[1 0.5 0 ];
f.LineWidth=1.6;
f.LineStyle='-';
legend({'\Omega_2 SOWFA','\Omega_2 simulated','\Omega_2 predicted','forecast frontier','\Omega_2 forecasted'},'Location','bestoutside','Orientation','horizontal')
legend('boxoff')

%dashed line to define frontier between data and prediction 

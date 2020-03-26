function [X1,X2,X3,X4,Y1,Y2,U1]=preprocessdmd2(beg, rotSpeed,time1,rotorAzimuth,nacelleYaw)

%% STATES


%states might be velocities or characteristics of the whole field.
%they are uploaded in the main script and chosen to serve as input to DMD
%function there

%% DETERMINISTIC STATES
%X1=omega1 ; X2=omega2; X3=omega1^2; X4=omega2^2;
X1=resample(rotSpeed(end-beg*10:1:end,1)',1,10);
X2=resample(rotSpeed(end-beg*10:1:end,2)',1,10);
% % X1_val=resample(detrend(rotSpeed_val(end-750*10:1:end,1)'),1,10);
% % X2_val=resample(detrend(rotSpeed_val(end-750*10:1:end,2)'),1,10);

X3=resample(rotSpeed(end-beg*10:1:end,1)',1,10).^2;
X4=resample(rotSpeed(end-beg*10:1:end,2)',1,10).^2;
% X3=X3./var(X3);
% X4=X4./var(X4);

%% OUTPUTS
%the rotor speeds of the two turbines are defined as outputs of the wind turbine system 
Y1=X1;
Y2=X2;

% % Y1_val=X1_val;
% % Y2_val=X2_val;

% Y1=resample(detrend(powerGenerator(end-750*10:1:end,1)'./1e6/rho),1,10);
% Y2=resample(detrend(powerGenerator(end-750*10:1:end,2)'./1e6/rho),1,10);

%% INPUTS: 
%BC: Multi-Blade Coordinate transformation 
%A directional thrust force  can be accomplished by implementinf MBC
%transformation, and decoupling/proejcting the blade loads in a non
%-rotating reference frame

% As a result, the measured out-of plane blade root bending moments M(t)
% --> [pitch{ij}(index,1);pitch{ij}(index,2);pitch{ij}(index,3)] are
% projected onto a non rotating reference frame --> PITCH 

%INPUTS: the inputs are the titl and yaw moments (fixed reference frame)

% Nturb=2;
% Offset=-8.4*2; %;-16;
% 
%     for index=1:1:length(time1)
%         %for each time instant INDEX get (for each turbine ij below) the 3
%         %out-of-plane blade root bending moments, corresponding to the
%         %three columns given a certain line
%        
%         for ij=1:1:Nturb
%             Azimuth=rotorAzimuth(index,ij);
%             PITCH=([1/3 1/3 1/3;
%                         2/3*cosd(Azimuth+Offset) 2/3*cosd(Azimuth+120+Offset) 2/3*cosd(Azimuth+240+Offset);  
%                         2/3*sind(Azimuth+Offset) 2/3*sind(Azimuth+120+Offset) 2/3*sind(Azimuth+240+Offset);])*...
%                     [pitch{ij}(index,1);pitch{ij}(index,2);pitch{ij}(index,3)];         
%                 
%             Pitch1(ij)=PITCH(1);
%             Pitch2(ij)=PITCH(2);
%             Pitch3(ij)=PITCH(3);
%             
%             %3 Matrixes containing the different bending moments where each
%             %line has the turbine number and each column the time instant
%             %INDEX
%             PPitch1(ij,index)=PITCH(1);
%             PPitch2(ij,index)=PITCH(2);
%             PPitch3(ij,index)=PITCH(3);
% 
% %             %Do the same for the validaion data
% %             Azimuth_val=rotorAzimuth_val(index,ij);
% %             PITCH=([1/3 1/3 1/3;
% %                         2/3*cosd(Azimuth_val+Offset) 2/3*cosd(Azimuth_val+120+Offset) 2/3*cosd(Azimuth_val+240+Offset);  
% %                         2/3*sind(Azimuth_val+Offset) 2/3*sind(Azimuth_val+120+Offset) 2/3*sind(Azimuth_val+240+Offset);])*[pitch_val{ij}(index,1);pitch_val{ij}(index,2);pitch_val{ij}(index,3)];
% %             
% %             Pitch1(ij)=PITCH(1);
% %             Pitch2(ij)=PITCH(2);
% %             Pitch3(ij)=PITCH(3);
% % 
% % 
% %             PPitch1_val(ij,index)=PITCH(1);
% %             PPitch2_val(ij,index)=PITCH(2);
% %             PPitch3_val(ij,index)=PITCH(3);
% 
%         end
%     end
    
i=1; 

%identification data
U1=resample(nacelleYaw(end-beg*10:1:end,1)',1,10);
% U1=U1./var(U1);
% U1=resample(detrend(PPitch2(i,end-750*10:1:end)),1,10); 
% U2=resample(detrend(PPitch3(i,end-750*10:1:end)),1,10);

%validation data
% % U1_val=resample(detrend(PPitch2_val(i,end-750*10:1:end)),1,10);
% % U2_val=resample(detrend(PPitch3_val(i,end-750*10:1:end)),1,10);
    

%Economic Model Predictive Control for Wake Steering: an Extended Dynamic Mode Decomposition Approach
%Master Thesis Dissertation
%Author: Nassir Rodrigues Cassamo
%Supervisors: Professor Jan Willem Van Wingerden and Professor João Sousa

%% RELEVAT INFORMATION
%olha aqui um grande testinho
%This script has the following goals:
% (1) Assess Simulation data, both qualitatively (animations) and
% quantitatively (graphics)
% (2) Derives a low dimensional model using Dynamica Mode Decomposition
% (variations included to take into account input-output data and other
% known states - deterministic states)
% (3) Validates the models with a set of validaiton data
% (4) Analyses the models, mainly the intrinsic dynamics  modes,
% damping, frequencies, energy) 
% (5) Reconstructs flow based on models an computes deviations from real 
% (6) Designs an Economic Model Predictive Control 

%This script requires the following functions to be added to MATLAB's path
% (1) Functions (LTI toolbox from DCSC, post processing tools from NRL
% (2) cbrewer: color scaling
% (3) altmany-export_fig: to export figures directly to specified directories

%This script requires data in the following fashion
% (1) 2 folders of name (2.1) steps_yaw and (2.2) steps_yaw_val, that
% contain the data non processed directly from CFD simulation from SOWFA
% (2) Data vectors with post processed information: the post processing is
% performed in the cluster and these vectors contain the resampled
% flowfield and resampled grid points. The grid points were resampled every
% Decimate (variable should be in these data vectors)
    % (2.1) U_data_complete_vec
    % (2.2) U_data_complete_vec_val

%% (0) INITIALISE
    %define type of simulation
     %pitch angle is varied and fixed frame loads will serve as inputs of model
    
    clc
    close all
    addpath('Functions')
    
    pitchmode=0;
    
    if pitchmode==0
        dirName={'steps_yaw'}; %directory for identification data: full results exported from SOWFA
        dirName_val={'steps_yaw_val'}; %directory for identification data: full resulsts exported from SOWFA
    elseif pitchmode==1
        dirName={'steps_theta_YT';}; %directory for identification data: full results exported from SOWFA
        dirName_val={'steps_theta_YT_val'}; %directory for identification data: full resulsts exported from SOWFA
    end
    videos=0;
    snapshots=0;

    % Turbine and flow characteristics to be used 
    rho=1.225; %air density in [kg m^-3]
    D=178; %Rotor Diameter used in simulations= 178 [m]
    % Simulation characterisitc (resampling)
    dt=2; %time sampling

%% (1) ASSESS DATA

    % IDENTIFICATION DATA
    visualisefirstresults(dirName,rho,0) %0: skip this; %1: see graphs
    % VALIDATION DATA
    visualisefirstresults(dirName_val,rho,0) %0: skip this; %1: see graphs

    if pitchmode==0
        filename='U_data_complete_vec_yaw.mat';
        load(filename) ;
        QQ_u=double(QQ_u);
        QQ_v=double(QQ_v);
    elseif pitchmode==1
        filename='U_data_complete_vec_pitch.mat';
        load(filename) ;
        QQ_u=full(QQ);
        QQ_u=double(QQ_u);
    end
 
    % Make movie from results (identification data, steady state (1500-3000)
     if videos==1
        %second input argument: dirctory to save movies
        [dirpathwake,dirpathveldef]=makeframes(D,'/Volumes/NASSIR/MATLAB/Movie_wake');
        [dirpathfinal]=makefinalframes(D,'/Volumes/NASSIR/MATLAB/Movie_combined');
        [dirpathpowerinsights]=makeframespowerinsights(D,rho,'/Volumes/NASSIR/MATLAB/power_insights');
        [dirpathvelfield]=plotvecfield(D,'/Volumes/NASSIR/MATLAB/power_velfield');
        [dirpathcuthubheightvec]=cuthubheightvec(D,'/Volumes/NASSIR/MATLAB/cuthubheight');
        
        %Movie making. Please specify Frames Per Second (FPS) inside
        makemoviewake('wake_deflection_yawcontrol_sowfa',dirpathwake);
        makemoviewake('wake_deflection_velocitydefifcitslice_yawcontrol_sowfa',dirpathveldef);
        makemoviewake('wake_deflection_combined',dirpathfinal);
        makemoviewake('wake_deflection_powerdelta',dirpathpowerinsights);
        makemoviewake('wake_deflection_velfield',dirpathvelfield);
        makemoviewake('wake_deflection_hubcut',dirpathcuthubheightvec);
        
        %Converts avi video to mp4. Also specify FPS inside
        movietomp4() 
     else
     end
     
    % Make figures to assess wake evolution (snapshots in sequential time
    % instants)
    if snapshots==1
        wake_vorticity_deflection
        wakevorticity_secondyaw
        hubheightcut
        cuthubtheightsecond
        vefield5D1
        velfield5D2
    else 
    end
 
%% (2) DYNAMIC MODE DECOMPOSITION 
    detrendingstates=1;
    begin=750;
    beg=750;

    % Read and process identification data
    [rotSpeed, nacelleYaw, time1,rotorAzimuth,pitch]=readdmdinformation(dirName); %read information from simulation
    [Inputs, Outputs, Deterministic]=preprocessdmd(beg, rotSpeed,time1,rotorAzimuth,nacelleYaw, pitchmode,pitch ); %preprocess information (resample and only relevant data)
    
    % Read and process validation data
    [rotSpeed_val, nacelleYaw_val, time1_val,rotorAzimuth_val,pitch_val]=readdmdinformation(dirName_val); %read information from simulation
    [Inputs_val, Outputs_val, Deterministic_val]=preprocessdmd(beg, rotSpeed_val,time1_val,rotorAzimuth_val,nacelleYaw_val,pitchmode,pitch_val); %preprocess information (resample and only relevant data)

    % Define states to be used for DMD
    states=QQ_u(:,(begin-beg)+1:end); % define states: first hypothesis 
    %states=[QQ_u(:,(begin-beg)+1:end);QQ_v(:,(begin-beg)+1:end);QQ_w(:,(begin-beg)+1:end)];

    if detrendingstates
        [states,meansteadystate,scalingfactor]=preprocessstates(states);
    else
    end
    
    r=120; %define truncation level for Singular Value Decomposition 
    maindir='/Volumes/NASSIR/MATLAB/'; %define directoty in user's computer to store all results
    [sys_red,FITje,U,S,V,method,X,X_p,Xd,dirdmd]=dynamicmodedecomposition(states,Inputs, Outputs, Deterministic,3,r,maindir); 

%% (3) DATA VALIDATION 
    % Validate Models from validation data set
    [FITje_val,dirdmd_val]=validatemodels(sys_red, Inputs_val,Outputs_val,r,strcat(dirdmd, '/val'));

    % Identification and Validation taks overview
    [modelVAF_val]=idvaloverview(FITje,FITje_val,dirdmd);

%% (4) DYNAMICAL ANALYSIS
    [f,LambdaDiag, P, phi,damping,b]=dynamicalanalysis(sys_red, U, S,V, dt,X_p,X, method,length(sys_red),0,D,9,Deterministic,r,dirdmd);
    visualisepodmodes(phi,f, P,x,y,z,Decimate,D,LambdaDiag,damping,method,Xd,dirdmd) 
    
%% (5) REBUILD FLOW FIELD AND ASSESS DEVIATIONS
    [statesrebuild]=rebuild(phi,b,LambdaDiag,r,X,Xd); %rebuild states with highest order model
  
    if detrendingstates
        for i=1:size(statesrebuild,2)
            statesrebuild(:,i)=statesrebuild(:,i)*scalingfactor+meansteadystate;
        end
    else
    end

    X=QQ_u(:,(begin-beg)+1:end-1);
    comparereconstruction(X, statesrebuild,D,dirdmd,x,y,z,Decimate,dirName)
    evauatemodelerror(X, statesrebuild,D,dirdmd,filename,dirName,x,y,z,Decimate)
    [statesfit]=evaluatetimevaryingerror(X,statesrebuild,dirdmd);
   
    save(strcat(dirdmd,'/RESULTS.mat'),'sys_red','FITje','FITje_val','X','statesrebuild','statesfit','Inputs','Outputs','Deterministic','Inputs_val','Outputs_val','Deterministic_val');
    close all
    
%% (6) ECONOMIC MODEL PREDICTIVE CONTROL DESIGN 
    [maxval,modeltouse]=max(FITje_val(2,:));
    [yf1, yf2]=evaluatepredictionpower(sys_red, modeltouse,Inputs_val, Outputs_val,dt, [1 350],400,200);



function [sys_red,FITje,U,S,V,method,X,X_p,Xd,dirdmd]=dynamicmodedecomposition(states, Inputs, Outputs, Deterministic ,method,r)

% This function aims to build a reduced order model from the states,
% input/output information and deterministic states gathered in the
% simulation and resampled

%% INPUT ARGUMENTS
%states: states to build DMD matrices
%inputs: [U1 U2 U3 ... Un];
%outputs: [Y1 Y3 .... Yn];
%deterministic: [Xd1 Xd2 Xd3 Xd4 ... Xdn];
%method: which method to use for DMD 
%r: truncation level (number of singular values to retain)


%% OUTPUTS ARGUMENTS
% sys_red: the state space systems, according to the number of singular
% values used
% FITje: the Fit of the model for the given data
% U, S, V: the matrices resulting from SVD to be used for DMD states reconstruction
%method: method is later used, as reconstruction depends on the method
%used, even though methodology is similar
% X, X_p: DMD matrices to be used later for reconstruction 

%% DMD - Dynamic Mode Decomposition

% x(k+1) ~ A x(k)
% X' ~ AX

% X = [   |  |        |  ]
%     [   x1 x2 ... xm-1 ]
%     [   |  |        |  ]

% X'= [   |  |        |  ]
%     [   x2 x3 ...  xm  ]
%     [   |  |        |  ]

%define necessary matrices for DMD
X      = states(:,1:end-1);
X_p    = states(:,2:end);

out      = Outputs(:,1:end-1);
out_p    = Outputs(:,2:end);

inp      =Inputs(:,1:end-1);
inp_p    =Inputs(:,2:end);

Xd     =Deterministic(:,1:end-1);
Xd_p   =Deterministic(:,2:end);

 
%% (1) DMDc - Dynamic Mode Decomposition with Control

if method==1 %dmdc algortihm
    
     dirdmd='/Volumes/NASSIR/MATLAB/DMDresults_DMDc';
    if ~exist(dirdmd,'dir') 
        mkdir(dirdmd);
    end
    
    %Goal of DMDc is to characterize the relationship between three
    %measurments: current measurment x(k), the future state x(k+1) and the
    %current control u(k). Relationship is approximated by the canonical
    %discrete linear dynamical system. 
    %Assumptions: all states are observable
    %are C will be identity.

    % x(k+1) ~  Ax(k) + Bu(k)
    % y(k)   ~  Cx(k) + Du(k)

    % X' ~ AX + BU 
    % Y  ~ CX + DU

    % X' ~ AX + BU ~ [A B][ X ] = G?
    %                     [ U ]

    %2 Perform SVD on the augmented data matrix SVD(?)=USV such that
    %G=X'VS^(-1)U*

    %3. Separate A and B by splitting left singular vectors into 2 separate
    %components

    % [A B] ~ [X'VS^(-1)U*1, X'VS^(-1)U*2 ]

    % U*1:
    % U*2:

    %4. Construct reduced order subspace from the output measurments X'. U
    %cannot be used to find low rank model of the dynamics and input matrixes
    %since it is defined for the input space which now includes both the state
    %measurments and the exogeneous inputs

    %SVD(X')=ÛS^V*

    % Ã  = Û*AÛ = Û*X'VS^(-1)U*1 Û
    % B~ = Û*B  = Û*X'VS^(-1)U*2
    
    Omega = [X;inp];
    
    [U, Sig, V]=svds(Omega,r);
    [Uf, Sigf, Vf]=svds(X_p,r);
    [Uo,So,Vo]=svds(X,r);
    FITje=zeros(2, r);
    OMEGA={};
    DAMPING={};

    for si=1:1:r
        
        Util=U(:,1:si);
        Sigtil=Sig(1:si,1:si);
        Vtil=V(:,1:si);
        
        Uhat=Uf(:,1:si);
        Sighat=Sigf(1:si,1:si);
        Vbar=Vf(:,1:si);
        
        n=size(X,1);
        q=size(Ups,1);
        U_1=Util(1:n,:);
        U_2=Util(n+q:n+q,:);
        
        approxA{si} = Uhat'*(X_p)*Vtil*inv(Sigtil)*U_1'*Uhat;
        approxB{si} = Uhat'*(X_p)*Vtil*inv(Sigtil)*U_2';
        approxC{si}=[Y1(:,1:end-1);Y2(:,1:end-1)]*pinv(([So(1:si,1:si)*Vo(:,1:si)'; U1(:,1:end-1)]));
        sys_red{si}=ss(approxA{si},approxB{si},approxC{si}(:,1:si),approxC{si}(:,si+1:end),2);
              
        [FITje,OMEGA,DAMPING]=evaluatemodel(sys_red,si,Inputs,Outputs,FITje,OMEGA,DAMPING,'identification');
        warning off
        export_fig(figure(1000+si),strcat(dirdmd,'/image',num2str(10000+si)),'-nocrop','-m2')
        warning on
    end
        
    close all
    VAFpermodes(FITje,r,{})
    warning off
    export_fig(figure(200),strcat(dirdmd,'/image',num2str(10000+si+1)),'-nocrop','-m2')
    warning on
           
  %% (2) ioDMD: Input Output Dynamic Mode Decomposition 
  
elseif method==2 %ioDMD
    
    %The goal of ioDMD is to capitalize on DMDc and extend it so a full
    %state space system may be obtained via usual subspace system
    %identification methods, estimating matrices A,B,C,D via lieast squares
    
    dirdmd='/Volumes/NASSIR/MATLAB/DMDresults_ioDMD';
    if ~exist(dirdmd,'dir') 
        mkdir(dirdmd);
    end
    
    dirdmdident='/Volumes/NASSIR/MATLAB/DMDresults_ioDMD/ident';
    if ~exist(dirdmdident,'dir') 
        mkdir(dirdmdident);
    end
    
    [U,S,V]=svds(X,r);
    FITje=zeros(2, r);
    OMEGA={};
    DAMPING={};

    for si=1:1:r
        
        Util=U(:,1:si);
        Sigtil=S(1:si,1:si);
        Vtil=V(:,1:si);
       
        all=[ Util'*X_p;out]*pinv([Sigtil*Vtil';inp]);
            
        A{si}=all(1:si,1:si);
        B{si}=all(1:si,si+1:end);
        C{si}=all(si+1:end, 1:si);
        D{si}=all(si+1:end, si+1:end);
        
        sys_red{si}=ss(A{si},B{si},C{si},D{si},2);
           
        [FITje,OMEGA,DAMPING]=evaluatemodel(sys_red,si,Inputs,Outputs,FITje,OMEGA,DAMPING,'identification');
        warning off
        export_fig(figure(1000+si),strcat(dirdmdident,'/image',num2str(10000+si)),'-nocrop','-m2')
        warning on
        close all
    end
        
    VAFpermodes(FITje,r,{})
    warning off
    export_fig(figure(200),strcat(dirdmdident,'/image',num2str(1000+length(sys_red)+1)),'-nocrop','-m2')
    warning on
    Xd={};
    
%% (3) extioDMD: Extended Input Output Dynamic Mode Decomposition
    
elseif method==3
    
    % The goal of extioDMD is to obtain a state space system, as ioDMD
    % performs, but etending the existing state to others which are not
    % necessarily related with the previous. This cpaitalizes on the
    % convergence of DMD results and the Koopman Operator, where there is
    % significant evidence that by incuding non linear observables (which
    % may be functions of the pre exisitng states, or not) DMD provides
    % better results.
    %These observables used to extend the current state space are referred
    %to as determinisitc states, as they are measurable and known for the
    %current scenario
    
    dirdmd='/Volumes/NASSIR/MATLAB/DMDresults_extioDMD';
        if ~exist(dirdmd,'dir') 
            mkdir(dirdmd);
        end
        
    dirdmdident='/Volumes/NASSIR/MATLAB/DMDresults_extioDMD/ident';
        if ~exist(dirdmdident,'dir') 
            mkdir(dirdmdident);
        end    
    
    [U,S,V]=svds(X,r);
    FITje=zeros(2, r);
    OMEGA={};
    DAMPING={};

    for si=1:1:r
        
        Util=U(:,1:si);
        Sigtil=S(1:si,1:si);
        Vtil=V(:,1:si);
       
        all=[ Xd_p; Util'*X_p;out]*pinv([Xd; Sigtil*Vtil';inp]);
            
        A{si}=all(1:size(Xd,1)+si,1:size(Xd,1)+si);
        B{si}=all(1:size(Xd,1)+si,size(Xd,1)+si+1:end);
        C{si}=all(size(Xd,1)+si+1:end, 1:size(Xd,1)+si);
        D{si}=all(size(Xd,1)+si+1:end, size(Xd,1)+si+1:end);
        
        sys_red{si}=ss(A{si},B{si},C{si},D{si},2);
        %same as before
           
        [FITje,OMEGA,DAMPING]=evaluatemodel(sys_red,si,Inputs,Outputs,FITje,OMEGA,DAMPING,'identification');
        warning off
        export_fig(figure(1000+si),strcat(dirdmdident,'/image',num2str(10000+si)),'-nocrop','-m2')
        warning on
        close all
    end
        
    VAFpermodes(FITje,r,{})
    warning off
    export_fig(figure(200),strcat(dirdmdident,'/image',num2str(1000+length(sys_red)+1)),'-nocrop','-m2')
    warning on
           
    
elseif method==4
    %% Professor Wingerden Least Square Solution for state space problem

    dirdmd='/Volumes/NASSIR/MATLAB/DMDresults_Wing';
    if ~exist(dirdmd,'dir') 
        mkdir(dirdmd);
    end
    
    % Take singular value decomposition of X with rank r
    % X ~ USV*
    % U belongs to set C with size nxr
    % S belongs to set C with sixe rxr
    % V belongs to set C iwth size mxr, and * denotes the conjugate transpose 
    % r is the rank of the reduced SVD approximation to X

    % left singular vectors U are POD modes
    % Columns of U are orthonormal, so U*U=I and V*V=I
   

    %Projection of full state space onto POD modes 
    %Make use of SVD decompositoin in this phase

    % U*X' ~ U*A(USV*) + U*BU 
    % Y  ~ C(USV*) + DU

    %A new state X^= U*X=U*USV*=SV*
    % Â  = U*A*U
    % B^ = U*B

    % U*X' ~ ÂSV* + B^U 

    % Matrixes Â and B^can now be found via least squares 

    % || U*X' - [ Â B^][ SV* ] ||
    %                  [  U  ]

    % [ Â B^ ] = U*X [ SV* ]* x  [ SV*VS    SV*U_*  ] ^-1
    %                [ U_   ]    [ U_SV     U_U_*   ]
  
    %%%---%%%---%%%
    
    %including deterministic states, the matrix problem formulation will
    %be:
    
    % [ I 0  ] [ Xd ] ~ [ I 0  ] A [ Xd   ] + [ I 0  ] BU
    % [ 0 U* ] [ X' ] ~ [ 0 U* ]   [ USV* ]   [ 0 U* ]
    
    % with 
    %
    % Â = [ I 0  ] A
    %     [ 0 U* ]
    
    % || [ I 0  ] [ Xd ]  -  Â[ Xd   ] - B^ U    || 
    % || [ 0 U* ] [ X' ]      [ USV* ]           || 
    
    % [ Â B^ ]
    
    %truncation/number of singular values for SVD decomposition
    
    Xd=[X1; X2; X3;X4];
    a=size(Xd);
    states=[states];
    [Uo,So,Vo]=svds(X,r);
    U=blkdiag(eye(size(Xd,1)),Uo); 
    QQ=[Xd; states];  
    X=QQ(:,1:end-1);
    X_p=QQ(:,2:end);
    S=blkdiag(eye(size(Xd,1)),So);
    V=[Xd(:,1:end-1)' Vo];
    FITje=zeros(2, r+a(1));
    OMEGA={};
    DAMPING={};
    for si=1:1:(r+a(1))
        
        Attt{si}=U(:,1:si)'*QQ(:,2:end)*[S(1:si,1:si)*V(:,1:si)';...
            [U1(:,1:end-1)]]'*inv([S(1:si,1:si)*V(:,1:si)'*V(:,1:si)*S(1:si,1:si) S(1:si,1:si)*V(:,1:si)'* [U1(:,1:end-1)]'; [U1(:,1:end-1)]*V(:,1:si)*S(1:si,1:si) [U1(:,1:end-1)]*[U1(:,1:end-1)]' ] );
        
        At{si}=Attt{si}(:,1:si);
        Bt{si}=Attt{si}(:,si+1:end);
        Ctt{si}=[Y1(:,1:end-1);Y2(:,1:end-1)]*pinv(([S(1:si,1:si)*V(:,1:si)'; U1(:,1:end-1)]));
        
        Ct{si} = Ctt{si}(:,1:si);
        Dt{si} = Ctt{si}(:,si+1:end);
        
        sys_red{si}=ss(At{si},Bt{si},Ct{si}, Dt{si},2);
        
        [FITje,OMEGA,DAMPING]=evaluatemodel(sys_red,si,Inputs,Outputs,FITje,OMEGA,DAMPING,'identification');
        
        warning off
        export_fig(figure(1000+si),strcat(dirdmd,'/image',num2str(10000+si)),'-nocrop','-m2')
        warning on
        close all
    end
        VAFpermodes(FITje,r,Xd)
        warning off
        export_fig(figure(200),strcat(dirdmd,'/image',num2str(10000+si+1)),'-nocrop','-m2')
        warning on
           
end    


    
end


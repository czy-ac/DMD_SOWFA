function [f,LambdaDiag, P, phi,damping,b]=dynamicalanalysis(sys_red, U, S,V, dt,X_p,X, method,mn,scaling,D,Uups,Deterministic,r,dirdmd)

%% DMD STANDARD MODAL ANALYSIS

if method==2 %ioDMD
    
    [W,Lambda]=eig(sys_red{mn}.A);
    LambdaDiag=diag(Lambda);
    omega=log(LambdaDiag)/dt/2/pi;
    f=abs(imag(omega))*D/Uups; %convertion from Hz to Strouhal number
    damping=-cos(atan2(imag(log(LambdaDiag)),real(log(LambdaDiag))));
    
    if scaling==0
        %phi=X_p*V*inv(S)*W;
        phi=U(:,1:mn)*W;
        b=phi\X(:,1);
        P=abs(b);
        
    elseif scaling==1 %common scaling by eigenvalues
        Ahat=(S(1:mn,1:mn)^(-1/2)) * sys_red{mn}.A * (S(1:mn,1:mn)^(1/2));
        [What,Dhat]=eig(Ahat);
        W_r=(S(1:mn,1:mn)^(1/2))*What;
        %Phi=U*W_r;
        phi=X_p*V(:,1:mn)/S(1:mn,1:mn)*W_r;
        P=(diag(phi'*phi));
    end

%% STATES HAVE BEEN EXTENDED TO ACCOUNT FOR NON LINEAR OBSERVABLES

elseif method==3 %extioDMD
    
    [W,Lambda]=eig(sys_red{mn}.A);
    LambdaDiag=diag(Lambda);
    omega=log(LambdaDiag)/dt/2/pi;
    f=abs(imag(omega))*D/Uups; %convertion from Hz to Strouhal number
    damping=-cos(atan2(imag(log(LambdaDiag)),real(log(LambdaDiag))));
    
    [sortedf,I]=sort(f);
    sorteddamping=damping(I);
    
    if scaling==0
        %phi=X_p*V*inv(S)*W;
        phi=blkdiag(eye(size(Deterministic,1),size(Deterministic,1)), U(:,1:mn))*W;
        b=phi\[Deterministic(:,1);X(:,1)];
        P=abs(b);
        
    elseif scaling==1
        Ahat=(S(1:mn,1:mn)^(-1/2)) * sys_red{mn}.A * (S(1:mn,1:mn)^(1/2));
        [What,Dhat]=eig(Ahat);
        W_r=(S(1:mn,1:mn)^(1/2))*What;
        %Phi=U*W_r;
        phi=X_p*V(:,1:mn)/S(1:mn,1:mn)*W_r;
        P=(diag(phi'*phi));  
    end

end

%% EIGENVALUE VISUALISATION 

% IN COMPLEX PLANE
figure460=figure('Position', [100 100 600 300]);
set(gcf,'color','w','Position', get(0, 'Screensize'));
subplot(1,2,1)
p=plot(LambdaDiag, 'o');
p.Color=[0.2 0.2 1];
p.MarkerSize=10;
p.MarkerEdgeColor=[0 0 0];
p.MarkerFaceColor=[.2 .2 1];
rectangle('Position', [-1 -1 2 2], 'Curvature', 1,'EdgeColor', 'k', 'LineStyle', '--');
axis(1.2*[-1 1 -1 1])
axis square
xlabel('Real axis \Re')
ylabel('Imaginary axis \Im')
title('Eigenvalue \lambda visualisation on the complex plane')
set(gca, 'FontSize', 14)
grid on
grid minor


% PER FREQUENCY OF OSCILLATION
subplot(1,2,2)
p2=plot(omega, 'o');
p2.Color=[0.2 0.2 1];
p2.MarkerSize=10;
p2.MarkerEdgeColor=[0 0 0];
p2.MarkerFaceColor=[.2 .2 1];
line([0 0], 0.3*[-1 1], 'Color', 'k','LineStyle', '--'); 
xlabel('Frequency \Omega [Hz]')
ylabel('Imaginary axis \Im')
title('Eigenvalue visualisation in frequencies per oscillation')
axis square
set(gca, 'FontSize', 14)
grid on
grid minor
export_fig(figure460,strcat(dirdmd,'/image','polescomplexplane'),'-nocrop','-m2');

%% POWER SPECTRUM 
 figure470=figure('Position', [100 100 600 300]);
 set(gcf,'color','w','Position', get(0, 'Screensize'));
 s=stem(f,P,'k','filled');
 grid on
 xlabel('Adimensionalised frequency St=f*D/U')
 ylabel('Power of mode || \phi ||')
 title('Dynamical Mode Decomposition Power Spectrum');
 set(gca, 'FontSize', 14)
 export_fig(figure470,strcat(dirdmd,'/image','dmdpowerspectrum'),'-nocrop','-m2');
 
 %% POWER ACCOUNTED FOR BY EACH MODE
 totalpower=sum(P);
 pfraction=P/totalpower;
 
 incpower=0;
 
 for i=1:length(P)
     incpower=incpower+pfraction(i);
     pp(i)=incpower;
 end
 
 figure(403)
 set(gcf,'color','w','Position', get(0, 'Screensize'));
 a=1:1:length(f);
 sid=scatter(a,pp,'o');
 hold on
 sid.MarkerFaceColor = [0 0.8 0.2];
 sid.MarkerEdgeColor = [0 0.8 0.2];
 powerplot=plot(a,pp,'LineWidth',1.5','color','red');
 powerplot.LineStyle='- -';
 powerplot.Color=[0 0.8 0.2];
 xlabel(' Number of singular values');
 ylabel('Total Energy maesured as || \phi || (%)')
 title('Energy accounted for by number of singular values used');
 set(gca, 'FontSize', 14)
 grid on
 grid minor
 hold off
 export_fig(figure(403),strcat(dirdmd,'/image','poweraccountedforwithsv'),'-nocrop','-m2');
 


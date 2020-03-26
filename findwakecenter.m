function [param, resnorm, residual, output,rebuildgaussian]=findwakecenter(i,D,x,y,z,Decimate,QQ_u)

%grid properties

Uups=9; %[m/s]

%% GRID PROPERTIES
figure
set(gcf,'color','w','Position', get(0, 'Screensize'));
[xx,yy,zz]=resamplegrid(x,y,z, Decimate);
yy=yy-500;
zz=zz-115;

X = length(xx);
Y = length(yy);
Z = length(zz);

k=70;
UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
Usecu=UmeanAbs_sh_u(:,k,:);
Usq=squeeze(Usecu);

[Ym_mesh , Zm_mesh]=meshgrid(yy/D,zz/D);
a=pcolor(Ym_mesh',Zm_mesh',(Usq-Uups)/Uups);
a.FaceAlpha=0.9;
shading interp
colormap(jet(4096))
view(2)
hold on
axis([ min(min(min(Ym_mesh))) max(max(max(Ym_mesh))) min(min(min(Zm_mesh))) max(max(max(Zm_mesh)))]);
xlabel('Distance [m]');
ylabel('Distance [m]');
pbaspect([1 1 1])
%[p]=plotrotor2D(0,D,0,115);
daspect([ 1 1 1])
set(gca,'fontsize', 12)

warning off
xxx=Ym_mesh;
yyy=Zm_mesh;
zzz=abs(Usq-Uups)/Uups;
[param,resnorm,residual,output] = gaussianwake(xxx,yyy,zzz);
warning on

s=scatter(param(4), param(5), 200, 'fill');
s.MarkerEdgeColor=[0.2 0.2 0.2];
s.MarkerFaceColor=[1 1 1];

%% reconstruct 3d wake and gaussian fit 
[xData, yData, zData] = prepareSurfaceData( xxx, yyy, zzz );
xyData = {xData,yData};

z = gaussian2Dapprox(param,xyData);
rebuildgaussian= reshape(z, [25 28]);











function [absVorState, absVelState,curlX]=getnewstates(QQ_u, QQ_v, QQ_w,x,y,z,Decimate)

[xx,yy,zz]=resamplegrid(x,y,z, Decimate);
[Xm_sh,Ym_sh,Zm_sh] = meshgrid(xx,yy,zz);
X = length(xx);
Y = length(yy);
Z = length(zz);

[l,c]=size(QQ_u);

for i=1:1:c
    
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    
    %curl_ circulation density of the fluid and provides mathematical
    %insights into fluid rotation based on fluid velocity field
    [CURLX, CURLY, CURLZ, CAV] = curl(Xm_sh,Ym_sh,Zm_sh,UmeanAbs_sh_u,UmeanAbs_sh_v,UmeanAbs_sh_w);
    
    absVor=sqrt(CURLX.^2+CURLY.^2+CURLZ.^2); %1
    absVel=sqrt(UmeanAbs_sh_u.^2+UmeanAbs_sh_v.^2+UmeanAbs_sh_w.^2); %2
    
    
    absVorState(:,i)=absVor(:);
    absVelState(:,i)=absVel(:);
    curlX(:,i)=CURLX(:);
    
end



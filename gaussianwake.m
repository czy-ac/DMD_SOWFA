function [param,resnorm,residual,output] = gaussianwake(xx,yy,zz)
% FMGAUSSFIT Create/alter optimization OPTIONS structure.
%   [fitresult,..., rr] = fmgaussfit(xx,yy,zz) uses ZZ for the surface 
%   height. XX and YY are vectors or matrices defining the x and y 
%   components of a surface. If XX and YY are vectors, length(XX) = n and 
%   length(YY) = m, where [m,n] = size(Z). In this case, the vertices of the
%   surface faces are (XX(j), YY(i), ZZ(i,j)) triples. To create XX and YY 
%   matrices for arbitrary domains, use the meshgrid function. FMGAUSSFIT
%   uses the lsqcurvefit tool, and the OPTIMZATION TOOLBOX. The initial
%   guess for the gaussian is places at the maxima in the ZZ plane. The fit
%   is restricted to be in the span of XX and YY.
%   See:
%       http://en.wikipedia.org/wiki/Gaussian_function
%          
%   Examples:
%     To fit a 2D gaussian:
%       [fitresult, zfit, fiterr, zerr, resnorm, rr] =
%       fmgaussfit(xx,yy,zz);
%   See also SURF, OMPTMSET, LSQCURVEFIT, NLPARCI, NLPREDCI.

%   Copyright 2013, Nathan Orloff.

%% Condition the data
[xData, yData, zData] = prepareSurfaceData( xx, yy, zz );
xyData = {xData,yData};

%% Set up the startpoint
[amp, ind] = max(zData); % amp is the amplitude.
xo = xData(ind); % guess that it is at the maximum
yo = yData(ind); % guess that it is at the maximum
sy = 10;
sx = 10;
xmax = max(xData);
ymax = max(yData);
xmin = min(xData);
ymin = min(yData);

%% Set up fittype and options.
Lower = [0,    0,   0,  xmin, ymin];
Upper = [Inf, Inf, Inf, xmax, ymax];
StartPoint = [amp, sx, sy, xo, yo];%[amp, sx, sy, xo, yo, zo];
% 
 tols = 1e-40;
 options = optimset('Algorithm','levenberg-marquardt',...
     'Display','off',...
     'MaxFunEvals',5e10,...
     'MaxIter',5e10,...
     'TolX',tols,...
     'TolFun',tols,...
     'TolCon',tols ,...
     'UseParallel','always');

%% perform the fitting
[param,resnorm,residual,exitflag, output, lambda,jacobian] = ...
    lsqcurvefit(@gaussian2D,StartPoint,xyData,zData,Lower,Upper);

end

function z = gaussian2D(par,xy)
% compute 2D gaussian
z=par(1)/(2*pi*par(2)*par(3))*...
    exp(-0.5*((xy{1}-par(4)).^2/(par(2).^2)+(xy{2}-par(5)).^2/(par(3).^2)));
end



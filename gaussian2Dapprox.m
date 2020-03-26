function z = gaussian2Dapprox(par,xy)
% compute 2D gaussian
z=par(1)/(2*pi*par(2)*par(3))*...
    exp(-0.5*((xy{1}-par(4)).^2/(par(2).^2)+(xy{2}-par(5)).^2/(par(3).^2)));
end



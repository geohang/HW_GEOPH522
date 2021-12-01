function rmse=model_variogram_error(h,V,c,a,type)
% variogram using the bounded linear and spherical models
% HPM 11/3/2020
% INPUT:
% h = lags for modeled estimates
% V = experimental variogram estimate
% c = variogram sill
% a = variogram range
% type = 'L' for linear and 'S' for spherical

Ix=find(h<=a); % finding lags less than a
switch type
    case 'L'
        Vm(Ix)=c*h(Ix)/a; % bounded linear
    case 'S'
        Vm(Ix)=c*(3*h(Ix)/(2*a)-0.5*(h(Ix)/a).^3); % spherical model
end

Ix2=h>a; % lags greater than range
Vm(Ix2)=c; % set equal to sill
V=V(:); Vm=Vm(:);
rmse=sqrt(mean((Vm-V).^2)); % root mean squared error


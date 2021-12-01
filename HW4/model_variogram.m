function V=model_variogram(h,c,a,type)
% variogram using the bounded linear and spherical models
% HPM 11/3/2020
% INPUT:
% h = lags for modeled estimates
% c = variogram sill
% a = variogram range
% type = 'L' for linear and 'S' for spherical

Ix=find(h<=a); % finding lags less than a
switch type
    case 'L'
        V(Ix)=c*h(Ix)/a; % bounded linear
    case 'S'
        V(Ix)=c*(3*h(Ix)/(2*a)-0.5*(h(Ix)/a).^3); % spherical model
end

Ix2=find(h>a); % lags greater than range
V(Ix2)=c; % set equal to sill

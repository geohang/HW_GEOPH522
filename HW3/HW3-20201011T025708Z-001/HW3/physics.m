function [rms]=physics(depth,vel,An)
% INPUTS:  x= independent variabe
%          y = dependent variable
%          An = parameters [A,n]
g=9.8; % [m/s^2]
rho=917; % [kg/m^3]
theta=10*pi/180; % convert slope angle to rad
um3=vel(1)-An(1).*(rho*g*sin(theta)).^An(2).*depth.^(An(2)+1); % Eq 6 in HW3
rms=sqrt(mean((um3-vel).^2)); % RMSE for each combo of n and A
end
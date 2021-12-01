function [rms]=physics1(depth,vel,u1,A)
% INPUTS:  x= independent variabe
%          y = dependent variable
%          An = parameters [A,n]
g=9.8; % [m/s^2]
rho=917; % [kg/m^3]
n=3;% use a fixed n
theta=10*pi/180; % convert slope angle to rad
um3=u1-A.*(rho*g*sin(theta)).^n.*depth.^(n+1); % Eq 6 in HW3
rms=sqrt(mean((um3-vel).^2)); % RMSE for each combo of n and A
end
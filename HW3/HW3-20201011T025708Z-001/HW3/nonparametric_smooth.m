function ymod = nonparametric_smooth(x,y,xmod,winsize)
% SNTX: ymod = nonparametric_smooth(x,y,xmod,winsize)
% this function smooths a 2-d dataset using a bisquare kernel
% INPUT: x = independent variable [n,1]
%        y = dependent variable [n,1]
%      xmod = locations of estimates [*,1]
%    winsize = size of the window (same units as x)
% OUTPUT ymod = nonparametric smoothed estimate
x=x(:); y=y(:); xmod=xmod(:);
ymod=zeros(size(xmod));
for i=1:length(xmod)
    dist=sqrt((x-xmod(i)).^2); % distance from each data point to the estimate location
    Ix=find(dist<winsize); % indicies to data within window
    Ix=Ix(isfinite(y(Ix)));% removing NaNs
    if isempty(Ix)
        ymod(i)=NaN; % use Nan if no data within window
    else
        w=15/16*(1-(dist(Ix)/winsize).^2).^2; % bisquare kernel
        ymod(i)=sum(w.*y(Ix))./sum(w); % unbiased estimate
    end
end
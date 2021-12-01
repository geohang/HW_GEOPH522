function [h,V] = semivariogram(x,y)
% simple 1D semivariogram function for equally spaced data
% INPUT:  
%  x = distance vector
%  y = measurement vector
% OUTPUT:
%  h = lag distance
%  V = semivariogram result
% SNTX: [h,V] = semivariogram(x,y)

% first define the lags
dx=mean(diff(x)); % average spacing
extent=(max(x)-min(x)); % extent
N=length(x); % number of data points
h=dx:dx:extent/2; % lags - only calculating to 1/2 extent to avoid bias
npairs=zeros(length(h),1); % preallocate number of pairs
V=zeros(length(h),1); % preallocate semivariance
for q=1:length(h) % loop over lags
    npairs(q)=N-q; % number of pairs at each lag
    Iu=1:(N-q); % index to heads
    Iv=(q+1):N; % index to tails
    V(q)=1/(2*npairs(q))*sum((y(Iu)-y(Iv)).^2); % semivariance
end 
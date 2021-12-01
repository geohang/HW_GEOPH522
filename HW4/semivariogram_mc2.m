function [h,V,npairs] = semivariogram_mc2(x,y,np)
% simple 1D semivariogram function for equally spaced data
% INPUT:  
%  x = distance vector
%  y = measurement vector
%  np = number of pairs of points to use 
% OUTPUT:
%  h = lag distance
%  V = semivariogram result
% SNTX: [h,V,npairs] = semivariogram_mc(x,y,np)

% first define the lags
dx=mean(diff(x)); % average spacing
extent=(max(x)-min(x)); % extent
N=length(x); % number of data points
h=dx:dx:extent/2; % lags - only calculating to 1/2 extent to avoid bias
npairs=zeros(length(h),1); % preallocate number of pairs
V=zeros(length(h),3); % preallocate semivariance
for q=1:length(h) % loop over lags
    npairs(q)=N-q; % number of pairs at each lag
    Iu=1:(N-q); % index to heads
    Iv=(q+1):N; % index to tails
    Vt=zeros(10,1);
    for m=1:100 % monte carlo for uncertainties
        I2=randsample(Iu,np); % random sampling of pairs
        Iut=Iu(I2);
        Ivt=Iv(I2);
        Vt(m)=1/(2*np)*sum((y(Iut)-y(Ivt)).^2); % semivariance
    end
    V(q,:)=quantile(Vt,[0.025 0.5 0.975]);
end 
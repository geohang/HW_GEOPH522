function [h,p,Vsave] = semivariogram_mc_model(x,y,np)
% simple 1D semivariogram function for equally spaced data
% INPUT:
%  x = distance vector
%  y = measurement vector
%  np = number of pairs of points to use
% OUTPUT:
%  h = lag distance
%  p = best parameters
%  Vsave = semivariogram result
% SNTX: [h,p,Vsave] = semivariogram_mc_model(x,y,np)

% first define the lags
dx=mean(diff(x)); % average spacing
extent=(max(x)-min(x)); % extent
N=length(x); % number of data points
h=dx:dx:extent/2; % lags - only calculating to 1/2 extent to avoid bias
npairs=zeros(length(h),1); % preallocate number of pairs
%V=zeros(length(h),3); % preallocate semivariance
% Vt=zeros(length(q),100);
for q=1:length(h) % loop over lags
    npairs(q)=N-q; % number of pairs at each lag
    Iu=1:(N-q); % index to heads
    Iv=(q+1):N; % index to tails
    
    for m=1:100 % monte carlo for uncertainties
        I2=randsample(Iu,np); % random sampling of pairs
        Iut=Iu(I2);
        Ivt=Iv(I2);
        Vt(q,m)=1/(2*np)*sum((y(Iut)-y(Ivt)).^2); % semivariance
    end
    
end

for m=1:100
    %for spherical model
    % first use brute force approach to find a suitable intital guess
    a=1:60;%brute force range for a
    c=0:60;%brute force range for c
    n=0:10;%brute force range for n
    
    for i=1:length(a)
        for j=1:length(c)
            for k=1:length(n)
                RMSE_3d1(i,j,k)=model_variogram_error_withnugget(h,Vt(:,m),c(j),a(i),n(k),'S');% calculate RMSE for the all a, c and n
            end
        end
    end
    %find the min RMSE and regarding index
    [min_val, position_min] = min(RMSE_3d1(:));
    [abestindex2,cbeatindex2,nbestindex2] = ind2sub(size(RMSE_3d1),position_min);
    
    abest2=a(abestindex2);%the optimum a in brute force method
    cbest2=c(cbeatindex2);%the optimum c in brute force method
    nbest2=n(nbestindex2);%the optimum n in brute force method
    
    
    
    fh=@(p)model_variogram_error_withnugget(h,Vt(:,m),p(1),p(2),p(3),'S');%function handle
    
    %The below method
%    [pbest fval ef]=fmincon(fh,[cbest2,abest2,nbest2],[],[],[],[],[min(h),min(Vt(:,m)),min(Vt(:,m))],[max(h) max(Vt(:,m)) max(Vt(:,m))]);
    [pbest,fval3]=fminsearch(fh,[cbest2,abest2,nbest2]);%gradient descent method to find best parameters
    p(m,1:3)=pbest;%store the best parameters
    V=model_variogram_withnugget(h,pbest(1),pbest(2),pbest(3),'S');
    Vsave(m,1:length(h))=V;%store the model variogram
% to test the plot
%     plot(h,Vt(:,m),'o')
%     hold on
%     plot(h,V)
%   %   text('c=21019900,a=1741897013,n=4.27')
%     set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
   
end

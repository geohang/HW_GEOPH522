% semivariograms, Lecture 20
D=load('DevonBdot.txt'); 
dist=D(:,1);
bdot=D(:,2);
figure(1);clf; 
plot(dist,bdot)
[h,V] = semivariogram(dist,bdot);
figure(2); clf
plot(h,V,'o')

%%
a=25; % range
c=45; % sill
Vm=model_variogram(h,c,a,'S')
figure(2); hold on
plot(h,Vm,'r','linewidth',2)

%% Brute Force approach to finding best variogram parameters
a=1:60;
c=30:60;
rmse=zeros(length(a),length(c));
for n=1:length(a)
    for m=1:length(c)
        rmse(n,m)=model_variogram_error(h,V,c(m),a(n),'L');
    end
end
figure(3);clf
imagesc(c,a,rmse);
colorbar

%% with gradient descent search
fh=@(p) model_variogram_error(h,V,p(1),p(2),'L');
[pbest,fval]=fminsearch(fh,[30,60])
figure(3); hold on
plot(pbest(1),pbest(2),'wo','markersize',8,'linewidth',2)







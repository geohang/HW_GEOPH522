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
[h,V,npairs] = semivariogram_mc(dist,bdot,100);
figure(3);clf
plot(h,V,'k')
figure(4);clf
plot(h,npairs)


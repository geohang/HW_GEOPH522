clear

%% load the data

dataset=load('DevonBdot.txt');

x=dataset(:,1);
y=dataset(:,2);

%% Q1

part1_x=x(1:100); %get the first section for x
part1_y=y(1:100); %get the first section for y
mean1(1)=mean(part1_y) % calculate the mean for part 1
std1(1)=std(part1_y) % calculate the standard deviation for part 1

part2_x=x(101:200);%get the second section for x
part2_y=y(101:200);%get the second section for y
mean1(2)=mean(part2_y) % calculate the mean for part 2
std1(2)=std(part2_y)% calculate the standard deviation for part 2

part3_x=x(201:300);%get the third section for x
part3_y=y(201:300);%get the third section for y
mean1(3)=mean(part3_y)% calculate the mean for part 3
std1(3)=std(part3_y)% calculate the standard deviation for part 3

part4_x=x(301:400);%get the fourth section for x
part4_y=y(301:400);%get the fourth section for y
mean1(4)=mean(part4_y)% calculate the mean for part 4
std1(4)=std(part4_y)% calculate the standard deviation for part 4

bins=30;
h=10;% window size for Kernel esimate

figure(1)
subplot(2,2,1)% for section 1
[centers] =plotRDH(part1_y,bins);%relative density histogram for the standard deviation value
hold on
[f] = myfun(part1_y,centers,h); % do the Kernel estimation 
plot(centers,f,'LineWidth',1.5)  % plot the Kernel estimation 
title('Section 1')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')

subplot(2,2,2)% for section 2
[centers] =plotRDH(part2_y,bins);%relative density histogram
hold on
[f] = myfun(part2_y,centers,h); % do the Kernel estimation 
plot(centers,f,'LineWidth',1.5)  % plot the Kernel estimation 
title('Section 2')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')

subplot(2,2,3)% for section 3
[centers] =plotRDH(part3_y,bins);%relative density histogram
hold on
[f] = myfun(part3_y,centers,h); % do the Kernel estimation 
plot(centers,f,'LineWidth',1.5)  % plot the Kernel estimation 
title('Section 3')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')

subplot(2,2,4)% for section 4
[centers] =plotRDH(part4_y,bins);%relative density histogram
hold on
[f] = myfun(part4_y,centers,h); % do the Kernel estimation 
plot(centers,f,'LineWidth',1.5)  % plot the Kernel estimation 
title('Section 4')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('Q1','-dpng')

%% Q3
[h,V] = semivariogram(x,y);% get the semivariance

var_data=var(y); %calculate the variance
cov_data=var_data-V; %calculate the covariance

auto_data=cov_data./var_data;%calculate the autocorrelation

figure(2); clf
plot(h,V,'o')%plot the semivariance
hold on
plot(h,cov_data,'linewidth',1)%plot the covariance
xlabel('lag')
legend('semivariance','covariance')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('Q3_1','-dpng')
figure(3)
plot(h,auto_data,'linewidth',1)%plot the autocorrelation
xlabel('lag')
ylabel('autocorrelation')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('Q3_2','-dpng')

%% Q4
np=10; % define the number of paris of points
[h1,V1,npairs1] = semivariogram_mc2(x,y,np);% semivariance for 10 paris

np=50; % define the number of paris of points
[h2,V2,npairs2] = semivariogram_mc2(x,y,np);% semivariance for 50 paris

np=100; % define the number of paris of points
[h3,V3,npairs3] = semivariogram_mc2(x,y,np);% semivariance for 100 paris

figure(4)
subplot(3,1,1)
plot(h1,V1(:,2),'k','linewidth',1)% plot the median of semivariance for 10 paris
hold on
plot(h1,V1(:,1),'b','linewidth',1)% plot the 95% uncertainty limits of semivariance for 10 paris
hold on
plot(h1,V1(:,3),'b','linewidth',1)% plot the 95% uncertainty limits of semivariance for 10 paris
legend('Median','95% uncertainty','95% uncertainty')
xlabel('lag')
ylabel('semivariance')
title('N_p=10')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')

subplot(3,1,2)
plot(h2,V2(:,2),'k','linewidth',1)% plot the median of semivariance for 50 paris
hold on
plot(h2,V2(:,1),'b','linewidth',1)% plot the 95% uncertainty limits of semivariance for 50 paris
hold on
plot(h2,V2(:,3),'b','linewidth',1)% plot the 95% uncertainty limits of semivariance for 50 paris
legend('Median','95% uncertainty','95% uncertainty')
xlabel('lag')
ylabel('semivariance')
title('N_p=50')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')

subplot(3,1,3)
plot(h3,V3(:,2),'k','linewidth',1)% plot the median of semivariance for 100 paris
hold on
plot(h3,V3(:,1),'b','linewidth',1)% plot the 95% uncertainty limits of semivariance for 100 paris
hold on
plot(h3,V3(:,3),'b','linewidth',1)% plot the 95% uncertainty limits of semivariance for 100 paris
legend('Median','95% uncertainty','95% uncertainty')
xlabel('lag')
ylabel('semivariance')
title('N_p=100')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('Q4','-dpng')

%% Q5

[h_part1,V_part1] = semivariogram(part1_x,part1_y);% semivariance for part 1

[h_part2,V_part2] = semivariogram(part2_x,part2_y);% semivariance for part 2

[h_part3,V_part3] = semivariogram(part3_x,part3_y);% semivariance for part 3

[h_part4,V_part4] = semivariogram(part4_x,part4_y);% semivariance for part 4

figure(5)

plot(h_part1,V_part1,'linewidth',1)%plot the semivariance for part 1
hold on
plot(h_part2,V_part2,'linewidth',1)%plot the semivariance for part 2
hold on
plot(h_part3,V_part3,'linewidth',1)%plot the semivariance for part 3
hold on
plot(h_part4,V_part4,'linewidth',1)%plot the semivariance for part 4
xlabel('lag')
ylabel('semivariance')
legend('Section 1','Section 2','Section 3','Section 4')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('Q5','-dpng')

%% Q6
% brute force approach

a=1:60; %brute force range for a
c=30:60; %brute force range for c
rmse=zeros(length(a),length(c)); %intial rmse

for n=1:length(a)
    for m=1:length(c)
        rmse(n,m)=model_variogram_error(h,V,c(m),a(n),'L'); % calculate RMSE for the all a and c
    end
end

figure(6);clf
imagesc(c,a,rmse);% image the rmse
xlabel('c')
ylabel('a')
colorbar
minvalue=min(min(rmse)) %find minimum RMSE
[abestindex,cbeatindex]=find(rmse==minvalue);%find index of minimum RMSE

abest=a(abestindex) %the optimum a in brute force method
cbest=c(cbeatindex) %the optimum c in brute force method
hold on
plot(cbest,abest,'o','MarkerSize',10,'MarkerFaceColor','k','linewidth',2); % plot the best point in image

print('Q6','-dpng')
%% Q7

fh=@(p)model_variogram_error(h,V,p(1),p(2),'L');%function handle
[pbest,fval]=fminsearch(fh,[30,60]);%gradient descent method to find best parameters


%% Q8

fh=@(p)model_variogram_error(h,V,p(1),p(2),'S');%function handle
[pbest1,fval1]=fminsearch(fh,[30,60]);%gradient descent method to find best parameters

%% Q9

%for linear model

% first use brute force approach to find a suitable intital guess
a=1:60;%brute force range for a
c=30:60;%brute force range for c
n=0:10;%brute force range for n
for i=1:length(a)
    for j=1:length(c)
        for k=1:length(n)
             RMSE_3d(i,j,k)=model_variogram_error_withnugget(h,V,c(j),a(i),n(k),'L');% calculate RMSE for the all a, c and n
        end
    end
end

%find the min RMSE and regarding index
[min_val, position_min] = min(RMSE_3d(:)); 
[abestindex1,cbeatindex1,nbestindex1] = ind2sub(size(RMSE_3d),position_min);

abest1=a(abestindex1)%the optimum a in brute force method
cbest1=c(cbeatindex1)%the optimum c in brute force method
nbest1=n(nbestindex1)%the optimum n in brute force method

%for linear model
fh=@(p)model_variogram_error_withnugget(h,V,p(1),p(2),p(3),'L');%function handle
[pbest_3d,fval2]=fminsearch(fh,[cbest1,abest1,nbest1]);%gradient descent method to find best parameters

%for spherical model

% first use brute force approach to find a suitable intital guess
a=1:60;%brute force range for a
c=30:60;%brute force range for c
n=0:10;%brute force range for n

for i=1:length(a)
    for j=1:length(c)
        for k=1:length(n)
             RMSE_3d1(i,j,k)=model_variogram_error_withnugget(h,V,c(j),a(i),n(k),'S');% calculate RMSE for the all a, c and n
        end
    end
end
%find the min RMSE and regarding index
[min_val, position_min] = min(RMSE_3d1(:)); 
[abestindex2,cbeatindex2,nbestindex2] = ind2sub(size(RMSE_3d1),position_min);

abest2=a(abestindex2)%the optimum a in brute force method
cbest2=c(cbeatindex2)%the optimum c in brute force method
nbest2=n(nbestindex2)%the optimum n in brute force method

%for spherical model
fh=@(p)model_variogram_error_withnugget(h,V,p(1),p(2),p(3),'S');%function handle
[pbest_3d1,fval3]=fminsearch(fh,[cbest2,abest2,nbest2]);%gradient descent method to find best parameters

%% Q10

V_linear=model_variogram(h,pbest(1),pbest(2),'L');%Linear model without nugget
V_linear_withnugget=model_variogram_withnugget(h,pbest_3d(1),pbest_3d(2),pbest_3d(3),'L');%Linear model with nugget

V_spherical=model_variogram(h,pbest1(1),pbest1(2),'S');%spherical model without nugget
V_spherical_withnugget=model_variogram_withnugget(h,pbest_3d1(1),pbest_3d1(2),pbest_3d1(3),'S');%spherical model with nugget

figure(7)
plot(h,V,'o') %plot the experimental variogram
hold on
plot(h,V_linear,'linewidth',1.5)%plot the Linear model without nugget
hold on
plot(h,V_spherical,'linewidth',1.5)%plot the spherical model without nugget
hold on
plot(h,V_linear_withnugget,'linewidth',1.5)%plot the Linear model with nugget
hold on
plot(h,V_spherical_withnugget,'linewidth',1.5)%plot the spherical model with nugget
legend('experimental variogram','Linear model without nugget','spherical model without nugget',...
    'Linear model with nugget','spherical model with nugget')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('Q10_1','-dpng')

figure(8)
subplot(2,2,1)
plot(h,V,'o') %plot the experimental variogram
hold on
plot(h,V_linear,'linewidth',1.5)%plot the Linear model without nugget
xlabel('lag')
ylabel('variogram')
title('Linear model without nugget')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
subplot(2,2,2)
plot(h,V,'o') %plot the experimental variogram
hold on
plot(h,V_spherical,'linewidth',1.5)%plot the spherical model without nugget
xlabel('lag')
ylabel('variogram')
title('spherical model without nugget')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
subplot(2,2,3)
plot(h,V,'o') %plot the experimental variogram
hold on
plot(h,V_linear_withnugget,'linewidth',1.5)%plot the Linear model with nugget
xlabel('lag')
ylabel('variogram')
title('Linear model with nugget')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
subplot(2,2,4)
plot(h,V,'o') %plot the experimental variogram
hold on
plot(h,V_spherical_withnugget,'linewidth',1.5)%plot the spherical model with nugget
xlabel('lag')
ylabel('variogram')
title('spherical model with nugget')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('Q10_2','-dpng')

%% Q11

%choose the spherical model with nugget which has smallest RMSE
np=50;%define the np
[h,psave,Vsave] = semivariogram_mc_model(x,y,np);%using new function to obtain the best parameters and variogram results

%% Q12

Vall=quantile(Vsave,[0.025 0.5 0.975]);% getmedian modeled variogram at each lag and 95% uncertainty limits
figure(9)
plot(h,Vall(2,:),'k','linewidth',1)% plot the median of semivariance 
hold on
plot(h,Vall(1,:),'b','linewidth',1)% plot the 95% uncertainty limits of semivariance
hold on
plot(h,Vall(3,:),'b','linewidth',1)% plot the 95% uncertainty limits of semivariance 
xlabel('lag')
ylabel('variogram')
legend('Median','95% uncertainty','95% uncertainty')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('Q11','-dpng')
%% Q13
bins=30;
h1=3;

figure(10)
subplot(1,3,1)
[centers] =plotRDH(psave(:,1),bins);%relative density histogram
hold on
[f] = myfun(psave(:,1),centers,h1); % do the Kernel estimation 
plot(centers,f,'LineWidth',1.5)  % plot the Kernel estimation 
title('sill')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')


subplot(1,3,2)
[centers] =plotRDH(psave(:,2),bins);%relative density histogram
hold on
[f] = myfun(psave(:,2),centers,h1); % do the Kernel estimation 
plot(centers,f,'LineWidth',1.5)  % plot the Kernel estimation 
title('range')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')


subplot(1,3,3)
[centers] =plotRDH(psave(:,3),bins);%relative density histogram
hold on
[f] = myfun(psave(:,3),centers,h1); % do the Kernel estimation 
plot(centers,f,'LineWidth',1.5)  % plot the Kernel estimation 
title('nugget')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')

print('Q13','-dpng')

%% Q14

%choose the spherical model with nugget which has smallest RMSE
np=150;%define the np
[h,psave,Vsave] = semivariogram_mc_model(x,y,np);%using new function to obtain the best parameters and variogram results

Vall=quantile(Vsave,[0.025 0.5 0.975]);% getmedian modeled variogram at each lag and 95% uncertainty limits
figure(11)
plot(h,Vall(2,:),'k','linewidth',1)% plot the median of semivariance 
hold on
plot(h,Vall(1,:),'b','linewidth',1)% plot the 95% uncertainty limits of semivariance
hold on
plot(h,Vall(3,:),'b','linewidth',1)% plot the 95% uncertainty limits of semivariance 
xlabel('lag')
ylabel('variogram')
legend('Median','95% uncertainty','95% uncertainty')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('Q14_1','-dpng')


bins=30;
h1=2;

figure(12)
subplot(1,3,1)
[centers] =plotRDH(psave(:,1),bins);%relative density histogram
hold on
[f] = myfun(psave(:,1),centers,h1); % do the Kernel estimation 
plot(centers,f,'LineWidth',1.5)  % plot the Kernel estimation 
title('sill')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')


subplot(1,3,2)
[centers] =plotRDH(psave(:,2),bins);%relative density histogram
hold on
[f] = myfun(psave(:,2),centers,h1); % do the Kernel estimation 
plot(centers,f,'LineWidth',1.5)  % plot the Kernel estimation 
title('range')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')


subplot(1,3,3)
[centers] =plotRDH(psave(:,3),bins);%relative density histogram
hold on
[f] = myfun(psave(:,3),centers,h1); % do the Kernel estimation 
plot(centers,f,'LineWidth',1.5)  % plot the Kernel estimation 
title('nugget')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')

print('Q14_2','-dpng')

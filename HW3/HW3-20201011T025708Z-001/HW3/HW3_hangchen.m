clear
clc

%% load the data
D = load('icevelocity.txt');%load the data
depth = D(:,1); %creates a column vector of depth
vel = D(:,2); %creates a column vector of velocites

%% Q1
degree=0:4;%initial the degree
vtest=zeros(length(vel),length(degree));%initial the vtest
for i=1:length(degree)
    P=polyfit(depth,vel,degree(i));% fit a line to the data;
    vtest(:,i)=polyval(P,depth);% evaluate at all depth
    rmse(i)=sqrt(mean((vtest(:,i)-vel).^2));% calculate the rmse
end

%codes for plot
figure;
clf
plot(depth,vel,'o') %plot the original data
hold on
plot(depth,vtest(:,1),'linewidth',1.5) %plot the data from models of degree 0
hold on
plot(depth,vtest(:,2),'linewidth',1.5) %plot the data from models of degree 1
hold on
plot(depth,vtest(:,3),'linewidth',1.5) %plot the data from models of degree 2
hold on
plot(depth,vtest(:,4),'linewidth',1.5) %plot the data from models of degree 3
hold on
plot(depth,vtest(:,5),'linewidth',1.5) %plot the data from models of degree 4
legend('Original',['RMSE(degree=0):',num2str(rmse(1))]...
    ,['RMSE(degree=1):',num2str(rmse(2))]...
    ,['RMSE(degree=2):',num2str(rmse(3))]...
    ,['RMSE(degree=3):',num2str(rmse(4))]...
    ,['RMSE(degree=4):',num2str(rmse(5))])%set the legend

xlabel('Depth (m)')% for the label of x axis
ylabel('Velocity (m/yr)')% for the label of y axis
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('Q1','-dpng')
%% Q2
%I use the getTrain function in this question for convenience
%The getTrain function is defind by class
degree=0:4;%initial the degree
pTrain=0.9;%define the pecent that's used in polynomial fit
nMC=1000; %times for Monte-Carlo
rmseCV=zeros(nMC,length(degree)); % initializing
for q=1:length(degree)
    for p=1:nMC
        [trainset, ~] = getTrainTest([depth vel],pTrain);%get 90% data
        ztrain=trainset(:,1); % depths for training
        vtrain=trainset(:,2); % velocity for training
        PP{q}(p,:)=polyfit(ztrain,vtrain,degree(q));% fit a line to the data;
        vm=polyval(PP{q}(p,:),ztrain); % evaluate at the train depths
        rmseCV(p,q)=sqrt(mean((vtrain-vm).^2));% calculate the RMSE
    end
end

%degree=0
mean_a0_0=mean(PP{1})%get the mean of A0 when degree=0
std_a0_0=std(PP{1})%get the standard deviation of A0 when degree=0
mean_rmse_0=mean(rmseCV(:,1))%get the mean of RMSE when degree=0
std_rmse_0=std(rmseCV(:,1))%get the standard deviation of RMSE when degree=0

%degree=1
mean_a1_1=mean(PP{2}(:,1))%get the mean of A1 when degree=1
std_a1_1=std(PP{2}(:,1))%get the standard deviation of A1 when degree=1
mean_a0_1=mean(PP{2}(:,2))%get the mean of A0 when degree=1
std_a0_1=std(PP{2}(:,2))%get the standard deviation of A0 when degree=1
mean_rmse_1=mean(rmseCV(:,2))%get the mean of RMSE when degree=1
std_rmse_1=std(rmseCV(:,2))%get the standard deviation of RMSE when degree=1

%degree=2
mean_a2_2=mean(PP{3}(:,1))%get the mean of A2 when degree=2
std_a2_2=std(PP{3}(:,1))%get the standard deviation of A2 when degree=2
mean_a1_2=mean(PP{3}(:,2))%get the mean of A1 when degree=2
std_a1_2=std(PP{3}(:,2))%get the standard deviation of A1 when degree=2
mean_a0_2=mean(PP{3}(:,3))%get the mean of A0 when degree=2
std_a0_2=std(PP{3}(:,3))%get the standard deviation of A0 when degree=2
mean_rmse_2=mean(rmseCV(:,3))%get the mean of RMSE when degree=2
std_rmse_2=std(rmseCV(:,3))%get the standard deviation of RMSE when degree=2

%degree=3
mean_a3_3=mean(PP{4}(:,1))%get the mean of A3 when degree=3
std_a3_3=std(PP{4}(:,1))%get the standard deviation of A3 when degree=3
mean_a2_3=mean(PP{4}(:,2))%get the mean of A2 when degree=3
std_a2_3=std(PP{4}(:,2))%get the standard deviation of A2 when degree=3
mean_a1_3=mean(PP{4}(:,3))%get the mean of A1 when degree=3
std_a1_3=std(PP{4}(:,3))%get the standard deviation of A1 when degree=3
mean_a0_3=mean(PP{4}(:,4))%get the mean of A0 when degree=3
std_a0_3=std(PP{4}(:,4))%get the standard deviation of A0 when degree=3
mean_rmse_3=mean(rmseCV(:,4))%get the mean of RMSE when degree=3
std_rmse_3=std(rmseCV(:,4))%get the standard deviation of RMSE when degree=3

%degree=4
mean_a4_4=mean(PP{5}(:,1))%get the mean of A3 when degree=4
std_a4_4=std(PP{5}(:,1))%get the standard deviation of A3 when degree=4
mean_a3_4=mean(PP{5}(:,2))%get the mean of A3 when degree=4
std_a3_4=std(PP{5}(:,2))%get the standard deviation of A3 when degree=4
mean_a2_4=mean(PP{5}(:,3))%get the mean of A2 when degree=4
std_a2_4=std(PP{5}(:,3))%get the standard deviation of A2 when degree=4
mean_a1_4=mean(PP{5}(:,4))%get the mean of A1 when degree=4
std_a1_4=std(PP{5}(:,4))%get the standard deviation of A1 when degree=4
mean_a0_4=mean(PP{5}(:,5))%get the mean of A0 when degree=4
std_a0_4=std(PP{5}(:,5))%get the standard deviation of A0 when degree=4
mean_rmse_4=mean(rmseCV(:,5))%get the mean of RMSE when degree=4
std_rmse_4=std(rmseCV(:,5))%get the standard deviation of RMSE when degree=4


%% Q3

%I use the getTrain function in this question for convenience
%The getTrain function is defind by class
%I also use the function defined plotRDH in HW2 to plot relative density hist

degree=0:4;%initial the degree
pTrain=0.9;%define the pecent that's used in polynomial fit
nMC=1000; %times for Monte-Carlo
rmseCV=zeros(nMC,length(degree)); % initializing
for q=1:length(degree)
    for p=1:nMC
        [trainset, testset] = getTrainTest([depth vel],pTrain);%get 90% data
        ztrain=trainset(:,1); % depths for training
        vtrain=trainset(:,2); % velocity for training
        ztest=testset(:,1);% depths for test
        vtest=testset(:,2);% velocity for test
        PPP{q}(p,:)=polyfit(ztrain,vtrain,degree(q));% fit a line to the data;
        vm=polyval(PPP{q}(p,:),ztest); % evaluate at the test depths
        rmseCV(p,q)=sqrt(mean((vtest-vm).^2));% calculate the RMSE
    end
end

bins=30;% set the bins
figure;

% plot the histgram
for i=1:5
    subplot(2,5,i)
    hist(rmseCV(:,i),bins) %degree from 0 to 4
    title(['degree=',num2str(i-1)])
    xlabel('Data value')% for the label of x axis
    ylabel('Counts')% for the label of y axis
    set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
end
% plot the relative density histgram
for i=1:5
    subplot(2,5,i+5)
    plotRDH(rmseCV(:,i),bins);%degree from 0 to 4
end
print('Q3','-dpng')
%% Q4

winsize=[3 10 50];% initial the windows' size
z0=0:2:180; %the points for the window average

for q=1:length(winsize)
    for p=1:length(z0)
        Ix=find(depth>(z0(p)-winsize(q)) & depth<(z0(p)+winsize(q))); % find values within window
        um(p,q)=mean(vel(Ix)); % take mean within window
    end
end
figure;
plot(depth,vel,'o') %plot the original data
hold on
plot(z0,um,'linewidth',2) %plot the moving window average data
legend('data','w=3','w=10','w=50')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('Q4','-dpng')
%% Q5
% I use the function nanparametric_smooth defined in class

winsize=[3 10 50];% initial the windows' size
z0=0:2:180;%the points for the window average

for q=1:length(winsize)
    ymod(q,:) = nonparametric_smooth(depth,vel,z0,winsize(q));%using a weighted moving window average
end

figure;
plot(depth,vel,'o') %plot the original data
hold on
plot(z0,ymod,'linewidth',2) %plot the moving window average data
legend('data','w=3','w=10','w=50')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('Q5','-dpng')

%% Q6
%first approach that I use the RMSE between the nonparametric models with
%physical models

pmodel=0.9;%define the pecent that's used in finding nonparametric models

g=9.8; % [m/s^2]
rho=917; % [kg/m^3]
theta=10*pi/180; % convert slope angle to rad
A=5e-18;%inital guess of A
n=3;%initial guess of n

um3=vel(1)-A.*(rho*g*sin(theta)).^n.*depth.^(n+1); % Eq 6 in HW3

winsize=[1:50];% initial the windows' size

[trainset, valiset] = getTrainTest([depth vel],pmodel);%get 90% data for model and 10% for validation

for q=1:length(winsize)
 
    
    ymod1= nonparametric_smooth(trainset(:,1),trainset(:,2),valiset(:,1),winsize(q));%using a weighted moving window average
    %, validation data for the location of esimation
       
    RMSE(q)=sqrt(mean((ymod1-valiset(:,2)).^2));% calculate the RMSE
end

[va,minloc]=min(RMSE) %output the minimum RMSE and locaiton

figure
plot(winsize,RMSE,'linewidth',2)%plot the winsize versus the RMSE
xlabel('Window size')
ylabel('RMSE')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('Q6_1','-dpng')

%second approach that I use the RMSE between the nonparametric models with
%measured data

% z0=0:2:180;%the points for the window average
% winsize=[1:50];% initial the windows' size
% 
% for q=1:length(winsize)
%     ymod1= nonparametric_smooth(depth,vel,z0,winsize(q));%using a weighted moving window average
%     RMSE(q)=sqrt(mean((ymod1-vel).^2));% calculate the RMSE
% end
% 
% 
% figure
% plot(winsize,RMSE)%plot the winsize versus the RMSE
% xlabel('Window size')
% ylabel('RMSE')
% set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
% print('Q6_2','-dpng')
%% Q7

% brute force approach
g=9.8; % [m/s^2]
rho=917; % [kg/m^3]
theta=10*pi/180; % convert slope angle to rad

n=2:0.01:4;%range of n
A=1e-18:0.1e-18:10e-18; %range of A

rms=zeros(length(A),length(n)); % initializing

for p=1:length(n)
    for q=1:length(A)
        um3=vel(1)-A(q).*(rho*g*sin(theta)).^n(p).*depth.^(n(p)+1); % Eq 6 in HW3
        rms(q,p)=sqrt(mean((um3-vel).^2)); % RMSE for each combo of n and A
    end
end

minvalue=min(min(rms)) %find minimum RMSE

[x y]=find(rms==minvalue)%find index of minimum RMSE

npt_n=n(y)%output the optimal n
npt_A=A(x)%output the optimal A

%% Q8

figure(7);clf
imagesc(n,A,rms,[0 10]); colorbar
xlabel('n')
ylabel('A')
hold on
plot(n(y),A(x),'o','MarkerSize',10,'MarkerFaceColor','k','linewidth',2);
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('Q8','-dpng')

%% Q9

fh=@(An)physics(depth,vel,An) % function handle, A can be tuned to data in v,z

A0=[A(x) n(y)];%give the initial value
[Abest,fval] = fminsearch(fh,A0) %use fminsearch to find optimal A and n

%% Q10

pTrain=0.9;%define the pecent that's used in polynomial fit
nMC=1000; %times for Monte-Carlo
rmseCV2=zeros(nMC,1); % initializing
u1=vel(1);

for p=1:nMC
    [trainset, ~] = getTrainTest([depth vel],pTrain);%get 90% data
    ztrain=trainset(:,1); % depths for training
    vtrain=trainset(:,2); % velocity for training
    fh=@(A)physics1(ztrain,vtrain,u1,A); % function handle, A can be tuned to data in v,z
    A0=Abest(1);% initial guess of A
    [Abest1,fval1] = fminsearch(fh,A0); %find the best parameters, and get the error
    Aall(p)=Abest1;% store the A
    rmseCV2(p)=fval1;% store the RMSE
end

bins=30;
figure
subplot(1,2,1)
plotRDH(Aall,bins);
title('A')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
subplot(1,2,2)
plotRDH(rmseCV2,bins);
title('RMSE')
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('Q10','-dpng')



%% Q11

figure(7)
hold on
%Plot the mean optimum values of A and its standard deviation with vertical errorbars
errorbar(n(y),mean(Aall),std(Aall),'+','linewidth',1)
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')

%% Q12

A_sample=mean(Aall)+randn(1,1000)*std(Aall);%generate 1000 simulated values
%with the mean and standard deviations of the A from my Monte-Carlo simulations.

R_sample=mean(rmseCV2)+randn(1,1000)*std(rmseCV2);%generate 1000 simulated values
%with the mean and standard deviations of the RMSE from my Monte-Carlo simulations.

%% Q13

[h,p] = kstest2(Aall,A_sample) %compare the actual distributions of A with
%simulated assuming a normal distribution

[h1,p1] = kstest2(rmseCV2,R_sample) %compare the actual distributions of RMSE with
%simulated assuming a normal distribution

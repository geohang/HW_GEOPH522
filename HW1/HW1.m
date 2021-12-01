clear all

% load the two datasets
D1=load('depths1.txt');% load the file depths1
D2=load('depths2.txt');% load the file depths2

%1.calculate and report the mean and standard deviation of two datasets
mean_depth1=mean(D1) %get the mean of depths1 and show it in command window
mean_depth2=mean(D2) %get the mean of depths2 and show it in command window

std_depth1=std(D1) %get the standard deviation of depths1 and show it in command window
std_depth2=std(D2) %get the standard deviation of depths1 and show it in command window

%2.calculate and report the median, mode, IQR, skewness, and kurtosis of two datasets
median_depth1=median(D1) %get the median of depths1 and show it in command window
median_depth2=median(D2) %get the median of depths2 and show it in command window

mode_depth1=mode(D1) %get the mode of depths1 and show it in command window
mode_depth2=mode(D2) %get the mode of depths2 and show it in command window

iqr_depth1=iqr(D1) %get the IQR of depths1 and show it in command window
iqr_depth2=iqr(D2) %get the IQR of depths2 and show it in command window

skewness_depth1=skewness(D1) %get the skewness of depths1 and show it in command window
skewness_depth2=skewness(D2) %get the skewness of depths2 and show it in command window

kurtosis_depth1=kurtosis(D1) %get the skewness of depths1 and show it in command window
kurtosis_depth2=kurtosis(D2) %get the skewness of depths2 and show it in command window

%3.boxplots for two datasets
figure;
boxplot(D1,'Notch','on','Labels','Depths1')% boxplot for depths1
axis([0 2 0 500]) %set axis limits
ylabel('Data value') % for the label of y axis
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('BOXPlot1','-dpng')

figure;
boxplot(D2,'Notch','on','Labels','Depths2')% boxplot for depths2
axis([0 2 0 500]) %set axis limits
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('BOXPlot2','-dpng')

% 4.histograms for two datasets
figure;
hist(D1,30)% histgram for depths1
xlabel('Data value')% for the label of x axis
ylabel('Counts')% for the label of y axis
axis([0 400 0 300]) %set axis limits
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('HISTPlot1','-dpng')
figure;
hist(D2,30)% histgram for depths2
axis([0 400 0 300]) %set axis limits
xlabel('Data value')% for the label of x axis
ylabel('Counts')% for the label of y axis
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('HISTPlot2','-dpng')


%5.boxplots and histograms for two datasets
figure
subplot(2,2,1)
boxplot(D1,'Notch','on','Labels','Depths1')% boxplot for depths1
axis([0 2 0 500]) %set axis limits
ylabel('Data value') % for the label of y axis
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
subplot(2,2,2)
boxplot(D2,'Notch','on','Labels','Depths2')% boxplot for depths2
axis([0 2 0 500]) %set axis limits
ylabel('Data value') % for the label of y axis
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
subplot(2,2,3)
hist(D1,30)% histgram for depths1
xlabel('Data value')% for the label of x axis
ylabel('Counts')% for the label of y axis
axis([0 400 0 300]) %set axis limits
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
subplot(2,2,4)
hist(D2,30)% histgram for depths2
axis([0 400 0 300]) %set axis limits
xlabel('Data value')% for the label of x axis
ylabel('Counts')% for the label of y axis
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('BOXHISTPlot','-dpng')

%6.boxplots and histograms for two datasets
figure
subplot(2,2,1)
boxplot(D1,'Notch','on','Labels','Depths1')% boxplot for depths1
axis([0 400 0 500]) %set axis limits
ylabel('Data value') % for the label of y axis
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
subplot(2,2,2)
boxplot(D2,'Notch','on','Labels','Depths2')% boxplot for depths2
axis([0 400 0 500]) %set axis limits
ylabel('Data value') % for the label of y axis
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
subplot(2,2,3)
hist(D1,30)% histgram for depths1
xlabel('Data value')% for the label of x axis
ylabel('Counts')% for the label of y axis
axis([0 400 0 500]) %set axis limits
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
subplot(2,2,4)
hist(D2,30)% histgram for depths2
axis([0 400 0 500]) %set axis limits
xlabel('Data value')% for the label of x axis
ylabel('Counts')% for the label of y axis
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('BOXHISTPlot1','-dpng')

%7.plot the relative density histogram of two datasets

figure
subplot(1,2,1)
[counts,centers] = hist(D1,30); %get the counts and centers of histogram of depths1
binWidth = centers(2)-centers(1);% calculate the width of each bin
bar(centers,counts/2000/binWidth) %plot the relative density histogram of depths1
hold on
axis([0 400 0 0.02]) %set axis limits
xlabel('Data value')% for the label of x axis
ylabel('PDF')% for the label of y axis
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
subplot(1,2,2)
[counts,centers] = hist(D2,30); %get the counts and centers of histogram of depths2
binWidth = centers(2)-centers(1);% calculate the width of each bin
bar(centers,counts/2000/binWidth) %plot the relative density histogram of depths2
axis([0 400 0 0.02]) %set axis limits
xlabel('Data value')% for the label of x axis
ylabel('PDF')% for the label of y axis
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('RDPlot','-dpng')

%8.plot Gaussian distribution with the relative density histogram
figure
subplot(1,2,1)
[counts,centers] = hist(D1,30); %get the counts and centers of histogram of depths1
binWidth = centers(2)-centers(1);% calculate the width of each bin
bar(centers,counts/2000/binWidth) %plot the relative density histogram of depths1
hold on

% plot Gaussian distribution with the relative density histogram in red
plot([0:400],mynormpdf([0:400],mean_depth1,std_depth1),'LineWidth',1,'Color','r')% line for dataset1
axis([0 400 0 0.02]) %set axis limits
xlabel('Data value')% for the label of x axis
ylabel('PDF')% for the label of y axis
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
subplot(1,2,2)
[counts,centers] = hist(D2,30); %get the counts and centers of histogram of depths2
binWidth = centers(2)-centers(1);% calculate the width of each bin
bar(centers,counts/2000/binWidth) %plot the relative density histogram of depths2
hold on

% plot Gaussian distribution with the relative density histogram in red
plot([0:400],mynormpdf([0:400],mean_depth2,std_depth2),'LineWidth',1,'Color','r')% line for dataset2
axis([0 400 0 0.02]) %set axis limits
xlabel('Data value')% for the label of x axis
ylabel('PDF')% for the label of y axis
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')
print('PDFPlot','-dpng')

%9.the probability of a new measurement at each site being within 20cm of the average value
%for depths1
pro_1_within20=length(D1(D1<(mean_depth1+20)&D1>(mean_depth1-20)))/2000 %get the probability
%for depths2
pro_2_within20=length(D2(D2<(mean_depth2+20)&D2>(mean_depth2-20)))/2000 %get the probability

%10.the probability of a new measurement at each site being at least 20cm larger than the average value
%for depths1
pro_1_out20=length(D1(D1>=(mean_depth1+20)))/2000 %get the probability
%for depths2
pro_2_out20=length(D2(D2>=(mean_depth2+20)))/2000 %get the probability

%11.the probability of a new measurement at each site being at least 20cm smaller than the average
%value
%for depths1
pro_1_smaller20=length(D1(D1<=(mean_depth1-20)))/2000 %get the probability
%for depths2
pro_2_smaller20=length(D2(D2<=(mean_depth2-20)))/2000 %get the probability

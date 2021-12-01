D = load('icevelocity.txt');
depth = D(:,1); %creates a column vector of depth
vel = D(:,2); %creates a column vector of velocites
figure(1);clf % open new figure
plot(depth,vel,'o','linewidth',2)
c = min(vel):1:max(vel); %creates variable c for FOR loop
for n = 1:length(c)
    chi_squared(n) = (1/length(vel))*sum((c(n)-vel).^2); % calculate cost function
end
figure(2);clf % new fig
plot(c,chi_squared,'linewidth',2) 
xlabel('C');
ylabel('\chi^2');
[chimin,Ix]=min(chi_squared); % find the minimum of the cost function
cbest=c(Ix); % best value of C
hold on; plot(c(Ix),chimin,'ro') % plot the best value
set(gca,'Linewidth',2,'fontsize',14) % make plot more readable
figure(1); hold on
cmodel=ones(length(vel),1)*cbest; % best model
plot(depth,cmodel,'r','linewidth',2)
set(gca,'Linewidth',2,'fontsize',14)
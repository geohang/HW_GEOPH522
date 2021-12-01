function [centers] = plotRDH(data,bins)
%function plotRDH
%input:
%     data: the data vector
%     bins: how many bins
%output:
%     centers: centers of histogram

[counts,centers] = hist(data,bins); %get the counts and centers of histogram 
binWidth = centers(2)-centers(1);% calculate the width of each bin
bar(centers,counts/length(data)/binWidth) %plot the relative density histogram
hold on
plot(centers,mynormpdf(centers,mean(data)...
    ,std(data)),'LineWidth',1.5) % plot the Gaussian curve for maximum value

xlabel('Data value')% for the label of x axis
ylabel('PDF')% for the label of y axis
set(gca,'LineWidth',1,'FontSize',14,'FontWeight','bold')


end


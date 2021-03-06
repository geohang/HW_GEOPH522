%%
D=load('elevations.txt');
D=D(:); %turns matrix into a column vector
x0=(nanmean(D)-nanstd(D)*10):nanstd(D)/10:(nanmean(D)+nanstd(D)*10);
N=hist(D,x0);
RDH=N/sum(N*nanstd(D)/10); % relative density histogram
figure(1);clf
bar(x0,RDH)

%%
h=20; % window size
[f] = myfun(D,x0,h);

% for n=1:length(x0)
%     dist=(D-x0(n)); % distance from x0 to all data values
%     Ix=find(abs(dist)<h); % finding all datapoints  within h of x0
%     w=15/16*(1-(dist(Ix)/h).^2).^2; % weights for all datapoints within h of x0
%     f(n)=nansum(w); % sum the weights
% end
% dx=nanstd(D)/10;
% f=1/sum(f*dx)*f; % normalized PDF so that it integrates to 1

%%
figure(1); hold on
plot(x0,f,'r','linewidth',2)




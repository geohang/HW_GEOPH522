function [f] = myfun(D,x0,h)
% it is used to generate kernel density estimate
% input D is the data, x0 is the central points and h is windows length



for n=1:length(x0)
    dist=(D-x0(n)); % distance from x0 to all data values
    Ix=find(abs(dist)<h); % finding all datapoints  within h of x0
    w=15/16*(1-(dist(Ix)/h).^2).^2; % weights for all datapoints within h of x0
    f(n)=nansum(w); % sum the weights
end
dx=nanstd(D)/10;
f=1/sum(f*dx)*f; % normalized PDF so that it integrates to 1
end


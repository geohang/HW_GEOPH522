 function f=mynormpdf(z,mu,sig)
 % f=mynormpdf(z,mu,sig)
 % INPUT: z = values to estimate PDF at
 % mu = mean value
 % sig = standard deviation
 % OUTPUT: f = probability density function for a normal/Gaussian distribution
 A=1/(sig*sqrt(2*pi)); % constant to ensure PDF integrates to 1 over +/? inf
 B=(z-mu).^2; % squared deviations from mean; numerator in exponential function
 C=2*sig.^2; % constant; denominator in exponential function
 f=A*exp(-B./C); % normal PDF
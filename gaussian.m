function y = gaussian(mu,sigma,n)
%function that calculates n random values y following the a Gaussian pdf
%  having a mean value mu and standard deviation sigma.
%
%   Inputs:
%    mu    = the mean of the Gaussian curve
%    sigma = the standard deviation
%    n     = length of the column vector

pi = 3.14159;

%create a column vector of uniformly distributed random numbers
randnow = -1+2.*rand([n,1]);       %random numbers between -1 and 1
y = mu + sqrt(2).*sigma.*erfinv(randnow);      %m_fines, a Gaussian


end
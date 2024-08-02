% X=ebisuzaki(x,nsim,value);
%
% This function creates 'nsim' random time series that have the same power
% spectrum as the original time series 'x' but random phases. The power
% spectrum of 'x' and 'X' is the same (see abs(fft(x)) and abs(fft(X));
% The 'X' time series has also the same variance as 'x'. (see Ebisuzaki W. (1997) ;
% A method to estimate the statistical significance of a correlation when 
% the data are serially correlated. J Climate, 10, 2147-2153).
%
% Input
% 'x' : vector of real number containing the time series to match
% 'nsim' : integer number giving the number of simulation
% 'value' : real number to initiate the random sequence (if > 0, the seeds
% is initiated to the value; otherwise, it is changed from the clock)
%
% Output
% 'X' : matrix of real number (the same number of rows as length of 'x' and
% 'nsim' columns.
%
% fixed by Hao Ye, 2012
%
% original version by
% Vincent MORON
% March 2002
function X = ebisuzaki(x,nsim,value)
if nargin<2
    nsim=1;
end
if nargin<3
    value=-1;
end
% 
% if(value>0)
%     rand('seed',value);
% else
%     rand('seed',sum(100*clock));
% end
if any(isnan(x))
    error('stats:ebisizaki','This procedure cannot handle NaNs');
end
    
n = length(x);
n2 = floor(n/2);
x = x(:);
a = fft(x); % discrete fourier transform of the original
amplitudes = abs(a); % power of the original spectrum
amplitudes(1) = 0;

for i = 1:nsim
   if (mod(n, 2) == 0) % series has even length
      thetas = 2*pi * rand(n2-1,1);
      angles = [0; thetas(:); 0; flipud(-thetas(:))];
      recf = amplitudes .* exp(1i*angles); % adjust the phases
      recf(n2,1) = sqrt(2) * recf(n2,1) * cos(rand()*2*pi); % adjust n/2 power
   else % series has odd length
      thetas = 2*pi * rand(n2, 1);
      angles = [0; thetas(:); flipud(-thetas(:))];
      recf = amplitudes .* exp(1i*angles); % adjust the phases
   end
   X(:,i) = real(ifft(recf));
end

% adjust variance of the surrogate time series to match the original
X = X * diag(std(x)./std(X));

% adjust mean of the surrogate time series to match the original
% X = X - repmat(mean(X), n, 1);
% X = X + mean(x);

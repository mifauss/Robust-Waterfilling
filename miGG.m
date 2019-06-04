function Mi = miGG( beta, h, P )

% INPUTS
   % beta- order of GG distribution 
   % P- power of the GG distribution
   % h - channel gain
   
% Output
   % Mi Mutual information 

% Note: 
%   To obtain more accurate result, the value of K and N can be
%   set larger


% Setting the integration range over y
K = 10000; 
Y = linspace(-7*abs(h*sqrt(P)),7*abs(h*sqrt(P)),K);   % integrating over 7 standard deviations


fy = zeros(1,K);


% number of random samples in MonteCarlo
N = 200000;


% Normalization Constants for GG
a = sqrt( gamma(1/beta)/gamma(3/beta)*P);

cbeta = beta/(2*a*gamma(1/beta));

for i = 1:K
    
    y = Y(i);
    
    Z = randn(1,N);
    
    fy(i) = mean(cbeta/abs(h)*exp(- abs(y-Z).^beta/(a*abs(h))^beta) );  % pdf of Y
    
end

Sum = trapz(Y,fy);

fy = fy/Sum;  % normalizing for it to be a proper density.


hY = -trapz(Y,fy.*log(fy));  %computing h(Y)

Mi = hY - 0.5*log(2*pi*exp(1));  % computing Mutual Information




end


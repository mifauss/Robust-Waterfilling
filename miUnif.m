function Mi = miUnif(h, P )

% INPUTS
   % P- power of the uniform distribution
   % h - channel gain
   
% Output
   % Mi Mutual information 

% Note: 
%   To obtain more accurate result, the value of K and N can be
%   set larger

    
% Setting the integration range over y
K = 8000; 
Y = linspace(-10,10,K);   % integrating over 10 standard deviations


fy = zeros(1,K);


% number of random samples in MonteCarlo
N = 8000; 


a = sqrt(3*P);


for i = 1:K
    
    y = Y(i);
    
    U = a*(rand(1,N)-0.5);
    
    fy(i) = mean( 1/sqrt(2*pi)*exp(- (y-h*U).^2/2) );  % pdf of Y
        
end

Sum = trapz(Y,fy);

fy = fy/Sum;  % normalizing for it to be a proper density.

hY = -trapz(Y,fy.*log(fy));  %computing h(Y)

Mi = hY - 0.5*log(2*pi*exp(1));  % computing Mutual Information



end


function [CallPrices, KK, k]=CarrMadden (funcHandle,parameters , S0,t,r,q)

% parameters chosen by Carr and madam for the pricing of the equity options
N = 4096;         
alpha = .75; %0.75;
grid_space = 0.25;  % space grid (iteration step) = eta in course notes
lambda = 2*pi/(N*grid_space);
b = lambda*N/2;
k = [-b:lambda:b-lambda]; % log-strike prices at which the value of the option are computed 

v = [0:grid_space:(N-1)*grid_space];
u = v-(alpha+1)*1i;

% S0=StockSpec(1);
% t=StockSpec(2);
% r=StockSpec(3);
% q=StockSpec(4);

% function handle 
    Characteristic=str2func(funcHandle);

rho = exp(-r*t)*Characteristic(parameters, S0,t,r,q,u)./(alpha^2+alpha-v.^2+1i*(2*alpha+1)*v);   
%correct the characteristic functions for VG 
 
    simpson_1 = (1/3);                      % to perform the Simpson integration 
    simpson = ((3 + (-1).^(2:1:N))/3);
    simpson_int = [simpson_1 simpson];    % dimension 1*N
    a = real(fft(rho.*exp( 1i*v*b)*grid_space.*simpson_int, N)); % Simpson rule


% CallPrices = exp(-alpha*k)/pi.*a;       
CallPrices = (1/pi)*exp(-alpha*k).*a ;      
% out=CallPrices;

KK = exp(k);
% out = spline(KK,CallPrices,K); % Cubic spline interpolation --> In y we will have the interpolated values of the CallPrices at the strike prices given in 'strikes'
% plot(KK(2500:1:3000),CallPrices(2500:1:3000))
% plot(KK,CallPrices)
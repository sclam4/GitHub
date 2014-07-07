function out = Hest_char(ModelParameters,S0,t,r,q,u)
%kappa = rate of reversion
%theta = long run variance >>>> eta
%sigma = Volatility of volatility >>>> lambda
%                                                                
%v0 = initial Variance
%rho = correlation between the two correlated brownian motions
%T = Time till maturity
%r = interest rate
%S0 = initial asset price

%at/in the money call
v0=ModelParameters(1);      %initial volatility
kappa=ModelParameters(2);   %mean reversion rate
eta=ModelParameters(3);     %long run variance
lambda=ModelParameters(4);  %vol of Vairance
rho=ModelParameters(5);     %correlation




% % % % alpha = 1.5;
% % % % N= 4096;
% % % % eta=0.25;
% % % % c=eta*N;
% % % % lambda=2*pi/c;
% % % % b =lambda*N/2;
% % % % u = [0:N-1]*eta; %u=eta*(j-1), p33 in world of VG
% % % % position = abs((log(strike) + b)/lambda + 1); %position of call - n in world of VG
%value in FFT
%matrix
% v = u - (alpha+1)*1i;
% zeta = -.5*(v.^2 +1i*v);
% gamma = kappa - rho*sigma*v*1i;
% PHI = sqrt(gamma.^2 - 2*sigma^2*zeta);
% A = 1i*v*(s0 + (r-q)*T);
% B = v0*((2*zeta.*(1-exp(-PHI.*T)))./(2*PHI - (PHI-gamma).*(1-exp(-PHI*T))));
% C = -kappa*theta/sigma^2*(2*log((2*PHI - (PHI-gamma).*(1-exp(-PHI*T)))./ (2*PHI)) + (PHI-gamma)*T);
% charFunc = exp(A + B + C); %in the notes world in VG this is noted by a phi
% ModifiedCharFunc = charFunc*exp(-r*T)./(alpha^2 + alpha - u.^2 + i*(2*alpha +1)*u); % in the notes world in VG this is noted by varrho (LaTeX)
% SimpsonW = 1/3*(3 + (-1).^[1:N] - [1, zeros(1,N-1)]); %slide 153 in advanced equity models
% FftFunc = exp(1i*b*u).*ModifiedCharFunc*eta.*SimpsonW; %World of VG - p 32 - alpha_j
% payoff = real(fft(FftFunc));
% CallValueM = exp(-log(strike)*alpha)*payoff/pi;
% format short;
% CallValue = CallValueM(round(position));


s0 = log(S0);

d=(     (rho*lambda.*u*1i-kappa).^2    -lambda^2.*(-1i*u-u.^2)     ).^0.5;
g=(kappa-rho*lambda.*u*1i-d) ./ (kappa-rho*lambda.*u*1i+d);

A= 1i.*u*(s0 + (r-q)*t);

B= eta*kappa*lambda^-2.*((kappa-rho*lambda.*u*1i-d)*t-2*log((1-g.*exp(-d*t))./(1-g)));

C=v0^2*lambda^-2.*(kappa-rho*lambda.*u*1i-d).*(1-exp(-d*t))./(1-g.*exp(-d*t));

out = exp(A + B + C);

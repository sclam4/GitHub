function y = BS_char (parameters, S0,t,r,q, u)
% S0=StockSpec(1);
% t=StockSpec(2);
% r=StockSpec(3);
% q=StockSpec(4);
sigma=parameters;
y = exp(1i*u*(log(S0)+(r-q-sigma^2/2)*t)).*exp(-sigma^2*t*u.^2/2);

function out=VG_char(parameters,S0,t,r,q,u) 
% paremeters= [C,G,M] ; 
% StockSpec=[S0,t,r,q]

C=parameters(1);
G=parameters(2);
M=parameters(3);
% S0=StockSpec(1);
% t=StockSpec(2);
% r=StockSpec(3);
% q=StockSpec(4);


omega= C*log((M-1)*(G+1)/(G*M));

out=exp(1i*u*(log(S0)+(r-q+omega)*t)).*(G*M./(G*M+(M-G)*1i.*u+u.^2)).^(C*t); %where is it????


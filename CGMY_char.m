function out=CGMY_char(para,S0,t,r,q,u)
C=para(1,1);
G=para(1,2);
M=para(1,3);
Y=para(1,4);
mu=r-q-C*gamma(-Y)*((M-1)^Y-M^Y+(G+1)^Y-G^Y);
out=exp(1i*u*(log(S0)+mu*t)).*exp(-t*C*gamma(-Y).*(M^Y-(M-1i.*u).^Y+G^Y-(G+1i*u).^Y));
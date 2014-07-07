function out=NIG_char(para, S0,t,r,q,u)
alpha=para(1);
beta=para(2);
delta=para(3);

mu=r-q+delta*(sqrt(alpha^2-(beta+1)^2)-sqrt(alpha^2-beta^2));
out=exp(1i*u*(log(S0)+mu*t)).*exp(-delta*t.*(sqrt(alpha^2-(beta+1i*u).^2)-sqrt(alpha^2-beta^2)));
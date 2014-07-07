% function [cd]=cdf2(FuncHandle, para, S0, t, r, q,x)

% parameter
para=OptimalParameters(1,:,2);
S0=169;
t=modifiedData(1,3,2);
r=modifiedData(1,4,2);
q=0;
FuncHandle=func2str(@VG_char)
%%
% find gridk, M, X0
chara=str2func(FuncHandle)
[call K k]=CarrMadden(FuncHandle, para, S0,t,r,q);
    X0_index=find(K>=0.001,1,'first')
    Xfinal_index=find(call>=0,1,'last')

    X0=k(X0_index)
    Xfinal=k(Xfinal_index)

    
    
Mtest=Xfinal_index-X0_index; 
if mod(Mtest,2)==0 
    M=Mtest/2; 
else M=(Mtest+1); end

theta=(Xfinal-X0)/(2*M); h=0.25; 
% 
% M=4096/2;
% theta=(Xfinal-X0)/(2*M)
% h=.25;

m=[-M:1:M];
v=m-.5;

char=chara(para,S0,t,r,q,v*h);


%
x=0;    
fm=exp(-1i*x.*v*h).*char./(v*pi);%.*simpson_int;

cd=1/2+1i/2*sum(fm)
function cd=cdfSumNegativea(FuncHandle, para, S0, t, r, q,kk)

% parameter
% para=OptimalParameters(1,:,2);
% S0=169;
% t=modifiedData(1,3,2);
% r=modifiedData(1,4,2);
% q=0;
% FuncHandle=func2str(@VG_char)
%%
if strncmpi(FuncHandle,'CGMY',4)
    C=para(1);G=para(2);M_=para(3);
    Y=.98;
    d_minus=-M_;
    d_plus=G;
%     k=exp(-t*C*gamma(-Y)*(M_^Y+G^Y)); 
    c=2*C*abs(gamma(-Y)*cos(pi*Y/2));
    nu=Y;
   
else c=4; 
end;



% find gridk, M, X0
chara=str2func(FuncHandle);
[call K k]=CarrMadden(FuncHandle, para, S0,t,r,q);
    X0_index=find(K>=0.001,1,'first');
    Xfinal_index=find(call>=0,1,'last');

    X0=k(X0_index);
    Xfinal=k(Xfinal_index);

    
%     
% Mtest=Xfinal_index-X0_index; 
% if mod(Mtest,2)==0 
%     M=Mtest/2; 
% else M=(Mtest+1); end

M=20;
delta=(Xfinal-X0)/(2*M); 
theta=delta/(2*pi);
a=(d_minus+d_plus)/2
d=2*min(-d_minus,d_plus);
h=(pi*d/(c*t))^(1/(1+nu))*M^(-nu/(1+nu));

m=[-M:1:M];
v=m*h+1i*a;

char=chara(para,S0,t,r,q,v);

    
fm=exp(-1i*X0*m*h).*char.*(1-(-1).^m.*exp(pi*a/h))./(v*pi/h);%.*simpson_int;

frft=frft2(fm,theta);

b=exp(a*(X0-delta.*kk)).*sum(frft);
cd=1+1i/2.*b;
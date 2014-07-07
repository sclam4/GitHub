function cd=cdfDIY(FuncHandle, para, S0, t, r, q)

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
   
elseif strncmpi(FuncHandle,'NIG',3)
    alpha=para(1);beta=para(2);deltaNIG=para(3);
     d_minus=beta-alpha;
      d_plus=beta+alpha;
%     k=exp(-t*delteNIG*sqrt(alpha^2-beta^2)); 
    c=deltaNIG;
    nu=1;
   
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

M=5;
gridk=(Xfinal-X0)/(2*M); 

a=(d_minus+d_plus)/2;
d=2*min(-d_minus,d_plus);
h=(pi*d/(c*t))^(1/(1+nu))*M^(-nu/(1+nu));
theta=gridk*h/(2*pi);

m=[-M:1:M];
v=m-.5;



kk=[0:1:2*M];


% find N and set up vector c and vector g*

l=fzero(@(x)2^x-2*(2*M+1)-1,10); 
N=2^(ceil(l))


cVector=zeros; mm=[0:1:2*M];mmInv=[2*M:-1:1];
cVector=[exp(1i*pi*mm.^2*theta) zeros(1,N-4*M-1)...
         exp(1i*pi*mmInv.^2*theta) ];
     

char=chara(para,S0,t,r,q,v*h);
fm=exp(-1i*X0*m*h).*char./(v); %.*simpson_int;
gVector=[exp(-1i*pi*mm.^2*theta).*fm zeros(1,N-2*M-1)];

% set up FMatrix and its inverse F(-1)=transpose(conj(F))/N (or F'/N)
FMatrix=zeros(N);
FMatrix(1,:)=ones(1,N);
w=exp(1i*2*pi/N);
for i=2:N for j=1:N, FMatrix(i,j)=w^((i-1)*(j-1)); end; end;
FInverse=FMatrix'/N;

%do two fourier transforms and one inverse fourier
A=FMatrix*transpose(cVector);
B=FMatrix*transpose(gVector);
ABInv=FInverse*(A.*B);


% do fft
FFT=fft(exp(-1i*pi*kk.^2*theta).*ABInv);

cd=1/2+1i/(2*pi)*exp(1i*X0*h/2)...
    *exp(1i*kk*gridk*h/2)...
    .*transpose(FFT);



% % ifourier(A.*B, )

% % frft=frft2(fm,theta)
% % cd=1/2+1i/2*exp(1i*gridk*h*kk/2).*frft'
% % 
% % 
% % % cd=Disfrft(fm,theta)
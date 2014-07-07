S0=169.18;
t=0.0795;
r=spline(yield(:,1), yield(:,2), t);
q=0;
initial=0.3;
data=modifiedData;
FunctionName=func2str(@BS_char);

Fit (0.3, FunctionName,data, S0,t,r,q)

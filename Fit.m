function [out, EstPrice]=Fit(ModelParameters, FunctionName, data, S0,q)

EstPrice=zeros;

DataStrike=data(:,1,:);
DataPrices=[data(1:20,2,1);data(:,2,2)];

% t= 0.07945 (29 days) , expiry date: 6th Sep 2013
t=data(1,3,1); r=data(1,4,1);
[CallPrices, KK]=CarrMadden(FunctionName, ModelParameters , S0, t,r,q);
for j=1:20
    EstPrice(j,1)=spline(KK,CallPrices,DataStrike(j,1,1));
end;
   
% % % t= 0.120548 (44 days) , expiry date: 21st Sep 2013
% % t=data(1,3,2); r=data(1,4,2);
% % [CallPrices, KK]=CarrMadden(FunctionName, ModelParameters , S0, t,r,q);
% % for j=1:95
% %     EstPrice(j,1,2)=spline(KK,CallPrices,DataStrike(j,1,2));
% % end;
% % 
% % EstPrices=[EstPrice(1:20,1,1);EstPrice(:,1,2)];
% % 


% t= 0.120548 (73 days) , expiry date: 19 Oct 2013
t=data(1,3,2); r=data(1,4,2);
[CallPrices, KK]=CarrMadden(FunctionName, ModelParameters , S0, t,r,q);
for j=1:56
    EstPrice(j,2)=spline(KK,CallPrices,DataStrike(j,1,2));
end;

EstPrices=[EstPrice(1:20,1);EstPrice(:,2)];


% least square function
y=(EstPrices-DataPrices).^2;
out=sum(y);





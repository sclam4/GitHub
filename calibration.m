% % OptimalParameters=zeros; %model*para*time
clc;clear;
S0=169.18; q=0;

options = optimset('PlotFcns',@optimplotfval , 'TolFun',1e-3,'Tolx',1e-3, 'MaxIter',1e6, 'MaxFunEvals', 1e6,'display','iter');



% %%% Calicrate BS
initial=.3;
makingdata();
data=modifiedData;
FunctionName=func2str(@BS_char);
tic;
[BSPara,fval]=fminsearchcon(@(para) Fit(para, FunctionName, data, S0,q),initial,[0],[],[],[],[], options);
time=toc;

OptimalParameters(1,:,1)=BSPara;
OptimalParameters(2,1,1)=time;
OptimalParameters(3,1,1)=fval;

% %%% Calibrate VG
initial=[1 5.87 5];
data=modifiedData;
FunctionName=func2str(@VG_char);
tic;
[VGPara,fval]=fminsearchcon(@(para) Fit(para, FunctionName, data, S0,q),initial,[0 0 0],[],[],[],[], options);
time=toc;
for j=1:3, OptimalParameters(1,j,2)=VGPara(j);end;
OptimalParameters(2,1,2)=time;
OptimalParameters(3,1,2)=fval;

% %%% Calibrate HEST
initial=[.0227 .2678 .1827 .3012 -.7520];
data=modifiedData;
FunctionName=func2str(@Hest_char);
tic;
[HestPara,fval]=fminsearchcon(@(para) Fit(para, FunctionName, data, S0,q),initial,[0 0 0 0 -inf],[],[],[],@(para) para(4)^2-2*para(2)*para(3) ,options);
time=toc;
for j=1:5, OptimalParameters(1,j,3)=HestPara(j);end;
OptimalParameters(2,1,3)=time;
OptimalParameters(3,1,3)=fval;

% %condition lamnda^2-2kappa*eta <0 and lower bound=[ 0 0 0 0 -inf]


%%% Calibrate CGMY
initial=[0.12943739856 8.074158819098 62.60142812967 0.99];
data=modifiedData;
FunctionName=func2str(@CGMY_char);
tic;
[CGMYpara,fval]=fminsearchcon(@(para) Fit(para, FunctionName, data, S0,q),initial,[0 0 0 0],[inf inf inf .999],[],[],[] ,options);
time=toc;
for j=1:4, OptimalParameters(1,j,4)=CGMYpara(j);end;
OptimalParameters(2,1,4)=time;
OptimalParameters(3,1,4)=fval;

%% Calibrate NIG
initial=[[33.10672393430 -20.92360355214 0.27404581997]];
data=modifiedData;
FunctionName=func2str(@NIG_char);
tic;
[NIGpara,fval]=fminsearchcon(@(para) Fit(para, FunctionName, data, S0,q),initial,[0 -inf 0],[],[],[],[],options);
time=toc;
for j=1:3, OptimalParameters(1,j,5)=CGMYpara(j);end;
OptimalParameters(2,1,5)=time;
OptimalParameters(3,1,5)=fval;


% % % % plot the graph with total error and calibration time
% % % % [CallPrice, KK]=CarrMadden ('VG_char', OptPara , S0, data(1,3,1),data(1,4,1),q);
% % % % [CallPrice2, KK2]=CarrMadden ('VG_char',OptPara , S0, data(1,3,2),data(1,4,2),q);

% % % % plot(DataStrike(:,1,1), DataPrice(:,1,1), 'bo', KK1(2700:1:2900), CallPrice1(2700:1:2900), 'rx',DataStrike(:,1,1), DataPrice(:,1,1), 'bo',KK2(2700:1:2900),CallPrice2(2700:1:2900),'g^')


DataStrike=[modifiedData(:,1,1) modifiedData(:,1,2)];
DataPrice=[modifiedData(:,2,1) modifiedData(:,2,2)];

subplot(3,2,1)

[x, EstPrice]=Fit(BSPara, 'BS_char', data, S0,q);
plot(DataStrike(1:20,1), DataPrice(1:20,1), 'bo', DataStrike(1:20,1), EstPrice(1:20,1), 'rx',...
    DataStrike(:,2), DataPrice(:,2), 'bo', DataStrike(:,2), EstPrice(:,2), 'rx');
title({['BS Model (least sq error: ',num2str(x),', time: ',num2str(OptimalParameters(2,1,1)),' sec)']})

subplot(3,2,2)

[x, EstPrice]=Fit(VGPara, 'VG_char', data, S0,q);
plot(DataStrike(1:20,1), DataPrice(1:20,1), 'bo', DataStrike(1:20,1), EstPrice(1:20,1), 'rx',...
    DataStrike(:,2), DataPrice(:,2), 'bo', DataStrike(:,2), EstPrice(:,2), 'rx');
title({['VG Model (least sq error: ',num2str(x),', time: ',num2str(OptimalParameters(2,1,2)),' sec)']})

subplot(3,2,3)

[x, EstPrice]=Fit(HestPara, 'Hest_char', data, S0,q);
plot(DataStrike(1:20,1), DataPrice(1:20,1), 'bo', DataStrike(1:20,1,1), EstPrice(1:20,1), 'rx',...
    DataStrike(:,2), DataPrice(:,2), 'bo', DataStrike(:,2), EstPrice(:,2), 'rx');
title({['HEST Model (least sq error: ',num2str(x),', time: ',num2str(OptimalParameters(2,1,3)),' sec)']})


subplot(3,2,4)

[x, EstPrice]=Fit(CGMYpara, 'CGMY_char', data, S0,q);
plot(DataStrike(1:20,1), DataPrice(1:20,1), 'bo', DataStrike(1:20,1,1), EstPrice(1:20,1), 'rx',...
    DataStrike(:,2), DataPrice(:,2), 'bo', DataStrike(:,2), EstPrice(:,2), 'rx');
title({['CGMY Model (least sq error: ',num2str(x),', time: ',num2str(OptimalParameters(2,1,4)),' sec)']})


subplot(3,2,5)

[x, EstPrice]=Fit(NIGpara, 'NIG_char', data, S0,q);
plot(DataStrike(1:20,1), DataPrice(1:20,1), 'bo', DataStrike(1:20,1,1), EstPrice(1:20,1), 'rx',...
    DataStrike(:,2), DataPrice(:,2), 'bo', DataStrike(:,2), EstPrice(:,2), 'rx');
title({['NIG Model (least sq error: ',num2str(x),', time: ',num2str(OptimalParameters(2,1,5)),' sec)']})




%% website: 
% http://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon/content/FMINSEARCHBND/demo/html/fminsearchcon_demo.html

%%
% a=1; b=2;
% x0=[2 3 4 5 6];
% lb=[1 2 3 4 5]; ub=[inf inf inf inf 8]; 
% [x,fval]=fminsearchcon(@(x) fit(x,a,b),x0,lb,ub,[],[],@(x) x(4)^2+x(3)^0.5-100)

% a=1; b=2;
% x0=[2 3 4 5 6];
% lb=[1 2 3 4 5]; ub=[inf inf inf inf 8]; 
% [x,fval]=fminsearchcon(@(x) fit(x,a,b),x0,lb,ub,[],[],@(x) x(4)^2+x(3)^0.5-100)

banana = @(x,a)a*(x(2)-x(1)^2)^2+(1-x(1))^2;
a=10; b=2;
x0=[4 100 ];
lb=[-5 2 ]; ub=[100 100]; 

options = optimset('PlotFcns',@optimplotfval , 'TolX', 1e-1,'TolFun',1e-3);
[x,fval]=fminsearchcon(@(x) banana(x,a),x0,[-5 -9],ub,[1 1; 1 2],[3000 5000],@(x) x(1)^3+x(2)^2-90000, options)
% myPlots.m
%
% f = 1/2 * (x1^2 + x2^2) + 1/2 * (x1+x2-1)^2 - (x1 + 2*x2);
% plot the changes of x1 x2 with different stopping criteria 
% |Xk-Xk-1| ranging from 10^-3 to 10^-8
%
%@Yibo Wang     yibow6@uci.edu
% 10-15-2015

clear all
clc

figure(1); %  plot average x1 vs stopping criteria  |Xk - Xk-1|
x = [10^-3, 10^-4, 10^-5, 10^-6, 10^-7, 10^-8];
y1 = [3.387403597834082e-01, 3.337069287968509e-01,3.332138396624784e-01, 3.333308396321729e-01, 3.333330974513671e-01, 3.333333285830687e-01];
loglog(x,y1,'LineWidth',2);
xlabel('stopping criteria |Xk - Xk-1|');
ylabel('average of x1, minimum position');

figure(2); %  plot variance and error bar of  x1 vs stopping criteria  |Xk - Xk-1|
x = [10^-3, 10^-4, 10^-5, 10^-6, 10^-7, 10^-8];
varx = [3.182144309671200e-03, 1.072640065855570e-05,8.768271110291447e-07,8.718212773732584e-11,8.217339983857404e-13,7.254495093571155e-15];
ebx = [5.641049822215010e-02, 3.275118419012616e-03,9.363904693177653e-04,9.337137020378668e-06,9.064954486293577e-07,8.517332383775542e-08];
%loglog(x,varx,'LineWidth',2);
L = varx-ebx;
U = varx+ebx;
errorbar(x,varx,L,U,'LineWidth',2);
xlabel('stopping criteria |Xk - Xk-1|');
ylabel('variance and errorbar of x, minimum position');
set(gca,'xscale','log');
set(gca,'yscale','log');

figure(3); % plot average f(x) vs stopping criteria 
x = [10^-3, 10^-4, 10^-5, 10^-6, 10^-7, 10^-8];
y3 = [-1.830233492757940e+00, -1.833315447262304e+00, -1.833332479956849e+00, -1.833333333245142e+00, -1.833333333332038e+00, -1.833333333333324e+00];
loglog(x,y3,'LineWidth',2);
xlabel('stopping criteria |Xk - Xk-1|');
ylabel('average of f(x), the optimum');

figure(4); % plot times of evaluation vs stopping criteria
x = [10^-3, 10^-4, 10^-5, 10^-6, 10^-7, 10^-8];
y4 = [21, 26, 30, 35, 39, 44];
loglog(x,y4,'LineWidth',2);
xlabel('stopping criteria |Xk - Xk-1|');
ylabel('average of times of evaluation');

figure(5); % average relative distance vs stopping criteria
x = [10^-3, 10^-4, 10^-5, 10^-6, 10^-7, 10^-8];
y5 = [3.099840575393769e-03, 1.788607102966466e-05, 8.533764850593428e-07, 8.819183161534739e-11, 1.295497042974603e-12, 8.437694987151190e-15];
loglog(x,y5,'LineWidth',2);
xlabel('stopping criteria |Xk - Xk-1|');
ylabel('average relative distance');


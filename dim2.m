% dim2.m
%
% this is a example of 2 dimension conjugate gradient method.
% number of samples is 50;
% the test function is 
% f = 1/2 * (x1^2 + x2^2) + 1/2 * (x1+x2-1)^2 - (x1 + 2*x2);
%
% mathematical optimium is f = -11/6  at x1=1/3, x2=4/3.
% the mathematical optimium is used in the mathematical analysis
%
% @author Yibo Wang         yibow6@uci.edu
% 10-15-2015

clear
clc
rng(10152015);
dim = 2;
sample = 50;
h = [1,2];
xall = rand(1,dim*sample)*1000 -500;

f = @(x) 0.5*(x(1)^2 + x(2)^2) + 0.5*(x(1)+x(2)-1)^2 - (h(1)*x(1)+h(2)*x(2));
g = @(x) [x(1)+x(1)+x(2)-1-h(1); x(2)+x(2)+x(1)-1-h(2)];

zero = [0;0];
for i=1:2:99
    tstart = tic;
    x = [xall(i);xall(i+1)]; %x0

    gOld = g(x); % g0
    dOld = -gOld; % d0
    k = 0;
    dis1=sqrt(x(1)^2+x(2)^2); dis2=0;
    while abs(dis1 - dis2)> 0.000001
        dis1 = dis2;
        
        gK = g(x);
        if k ==0
            betaOld = 1;
        else
            betaOld = (gK' * gK)/ (gOld' * gOld);
        end
        dK = -gK + betaOld.* dOld;
        alphaK = (gK' * gK) / (dK' * (g(dK) - g(zero)));
        x = x + alphaK .* dK;
        
        dis2 = sqrt(x(1)^2+x(2)^2);
        k=k+1;
        gOld = gK;
        dOld = dK;
    end
    t(i) = toc(tstart);
    % record all the results
  
    for ii=1:2
        xresult(i+ii-1) = x(ii);        
    end
    toe(i) = k;
    fresult(i) = f(x);
    
    
end

%
% mathematical analysis 
%
sumX1=0; sumX2=0; sumF=0; sumToe=0; sumT=0; sumDis=0;
optimum = -11/6;
for i=1:2:99
    sumX1 = sumX1+xresult(i);
    sumX2 = sumX2+xresult(i+1);
    sumF = sumF+fresult(i);
    sumToe = sumToe+toe(i);
    sumT = sumT+t(i);
    sumDis = sumDis + abs(fresult(i) - optimum);
end
avrX1 = sumX1/sample;
avrX2 = sumX2/sample;
avrF = sumF/sample;
avrToe = sumToe/sample;
avrT = sumT/sample;
avrDis = sumDis/sample;

varX1=0; varX2=0; varF=0; varToe=0; varT=0; varDis=0;
for i=1:2:99
    varX1 = varX1 + (xresult(i)-avrX1)^2;
    varX2 = varX2 + (xresult(i+1)-avrX2)^2;
    varF = varF +(fresult(i)-avrF)^2;
    varToe = varToe + (toe(i)-avrToe)^2;
    varT = varT + (t(i)-avrT)^2;
    varDis = varDis + (abs(fresult(i)-optimum)-avrDis)^2;
end
varX1 = varX1/(sample-1);
varX2 = varX2/(sample-1);
varF = varF/(sample-1);
varToe = varToe/(sample-1);
varT = varT/(sample-1);
varDis = varDis/(sample-1);

ebX1 = sqrt(varX1);
ebX2 = sqrt(varX2);
ebF = sqrt(varF);
ebToe = sqrt(varToe);
ebT = sqrt(varT);
ebDis = sqrt(varDis);


% save the results to file 2d.txt
fid = fopen('2d.txt','wt'); 
fprintf(fid,'the average of x1 is %.15e\n',avrX1);  
fprintf(fid,'the variance of x1 is %.15e\n',varX1);  
fprintf(fid,'the error bar of x1 is %.15e\n',ebX1);

fprintf(fid,'the average of x2 is %.15e\n',avrX2);
fprintf(fid,'the variance of x2 is %.15e\n',varX2); 
fprintf(fid,'the error bar of x2 is %.15e\n',ebX2);  

fprintf(fid,'the average of f is %.15e\n',avrF);
fprintf(fid,'the variance of f is %.15e\n',varF);
fprintf(fid,'the error bar of f is %.15e\n',ebF);

fprintf(fid,'the average of times of evaluation is %.15e\n',avrToe);
fprintf(fid,'the variance of times of evaluation is %.15e\n',varToe);
fprintf(fid,'the error bar of times of evaluation is %.15e\n',ebToe);

fprintf(fid,'the average of running time is %.15e\n',avrT);
fprintf(fid,'the variance of running time is %.15e\n',varT);
fprintf(fid,'the error bar of running time is %.15e\n',ebT);

fprintf(fid,'the average of relative distance is %.15e\n',avrDis);
fprintf(fid,'the variance of relative distance is %.15e\n',varDis);
fprintf(fid,'the error bar of relative distance is %.15e\n',ebDis);

fclose(fid);

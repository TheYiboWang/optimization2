% dim10.m 
%
% this is a example of 10 dimension conjugate gradient method.
% number of samples is 50, samples are chosen from (-500,500);
% 
% the test function is
% f = 1/2 * sum(Xi^2) + 1/2 * (sum(Xi)-1)^2 - (x1 + 2*x2 + ... +10*x10);
%
% mathematical optimium is f = -60  at (-4,-3,-2,-1,0,1,2,3,4,5)
% the mathematical optimium is used in the mathematical analysis
%
% @author Yibo Wang         yibow6@uci.edu
% 10-15-2015

clear
clc
rng(10152015);
dim = 10;
sample = 50;
h = [1,2,3,4,5,6,7,8,9,10];
xall = rand(1,dim*sample) * 1000 - 500; %

f = @(x) 0.5*(x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2 + x(5)^2 + x(6)^2 + x(7)^2 + x(8)^2 + x(9)^2 + x(10)^2) + 0.5*(x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+x(10)-1)^2 - (h(1)*x(1)+h(2)*x(2)+h(3)*x(3)+h(4)*x(4)+h(5)*x(5)+h(6)*x(6)+h(7)*x(7)+h(8)*x(8)+h(9)*x(9)+h(10)*x(10));
g = @(x) [x(1)+x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+x(10)-1-h(1);
    x(2)+x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+x(10)-1-h(2);
    x(3)+x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+x(10)-1-h(3);
    x(4)+x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+x(10)-1-h(4);
    x(5)+x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+x(10)-1-h(5);
    x(6)+x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+x(10)-1-h(6);
    x(7)+x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+x(10)-1-h(7);
    x(8)+x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+x(10)-1-h(8);
    x(9)+x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+x(10)-1-h(9);
    x(10)+x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+x(10)-1-h(10)];

zero = zeros(dim,1); x = zeros(dim,1);
for i=1:dim:dim*sample-1
    tstart = tic;
    for ii=1:dim
        x(ii,1) = xall(i+ii-1); %x0
    end
    
    fTmp = 0;
    gOld = g(x); % g0
    dOld = -gOld; % d0
    k = 0;
    while abs(fTmp - f(x))> 0.000001
        fTmp = f(x);
        
        gK = g(x);
        if k == 0
            betaOld = 1;
        else
            betaOld = (gK' * gK)/ (gOld' * gOld);
        end
        dK = -gK + betaOld.* dOld;
        alphaK = (gK' * gK) / (dK' * (g(dK) - g(zero)));
        x = x + alphaK .* dK;
        
        k=k+1;
        gOld = gK;
        dOld = dK;
    end
    t(i) = toc(tstart);
    % record all the results
    
    for ii=1:dim
        xresult(i+ii-1) = x(ii);
    end
    toe(i) = k;
    fresult(i) = f(x);
    
    
end

%
% mathematical analysis
%
sumX=zeros(1,dim); avrX=zeros(1,dim);
sumF=0; sumToe=0; sumT=0; sumDis=0; optimum = -60;
for i=1:dim:dim*sample-1
    for ii=1:dim
        sumX(ii) = sumX(ii)+xresult(i+ii-1);
    end
    sumF = sumF+fresult(i);
    sumToe = sumToe+toe(i);
    sumT = sumT+t(i);
    sumDis = sumDis + abs(fresult(i) - optimum);
end
for ii=1:dim
    avrX(ii) = sumX(ii)/sample;
end
avrF = sumF/sample;
avrToe = sumToe/sample;
avrT = sumT/sample;
avrDis = sumDis/sample;

varX=zeros(1,10); ebX=zeros(1,10); varF=0; varToe=0; varT=0; varDis=0;
for i=1:dim:dim*sample-1
    for ii=1:dim
        varX(ii) = varX(ii) + (xresult(ii+i-1)-avrX(ii))^2;
    end
    varF = varF +(fresult(i)-avrF)^2;
    varToe = varToe + (toe(i)-avrToe)^2;
    varT = varT + (t(i)-avrT)^2;
    varDis = varDis + (abs(fresult(i)-optimum)-avrDis)^2;
end
for ii=1:dim
    varX(ii) = varX(ii)/(sample-1);
end
varF = varF/(sample-1);
varToe = varToe/(sample-1);
varT = varT/(sample-1);
varDis = varDis/(sample-1);

for ii=1:dim
    ebX(ii) = sqrt(varX(ii));
end
ebF = sqrt(varF);
ebToe = sqrt(varToe);
ebT = sqrt(varT);
ebDis = sqrt(varDis);


% save the results to file 10d.txt
fid = fopen('10d.txt','wt');
fprintf(fid,'the average of x is %.15e\n',avrX);
fprintf(fid,'the variance of x is %.15e\n',varX);
fprintf(fid,'the error bar of x is %.15e\n',ebX);

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

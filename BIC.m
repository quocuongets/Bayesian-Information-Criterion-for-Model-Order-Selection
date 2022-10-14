clear all;
close all;
clc;

load signals.mat;
y = Ay_mod;
x = Fy_mod;
N=size(y,1);
fs = 512;
outputs=size(y,2);

TT=0:1/fs:(length(y)-1)/fs;  % time (second)

% Plot data

figure (1)
subplot(2,1,1)
plot(TT,x)
xlabel('Time (s)','Fontsize',1.5)
ylabel('Force (N)','Fontsize',1.5)
xlim([0 TT(end)])
subplot(2,1,2)
plot(TT,y)
xlabel('Time (s)','Fontsize',1.5)
ylabel('Force (N)','Fontsize',1.5)
xlim([0 TT(end)])

% Histogram and Normal Probability Plot

figure (2)
subplot(1,2,1),histfit(y,round(sqrt(N)))
xlabel('Residuals')
ylabel('Probability')
title('Historgram of Acceleration output')

subplot(1,2,2), normplot(y)
xlabel('Data')
ylabel('Probability')
title('Normal probability plot')

% Plot ACF and PACF
figure (3)
subplot(1,2,1)
autocorr(y,'NumLags',100)
subplot(1,2,2)
parcorr(x,'NumLags',100)

% ARX model

DATA=iddata(y,x,1/fs);
    
minar=2;
maxar=100;
tic
        for order=minar:maxar
            models{order}=arx(DATA,[order order 0]);
        end
toc
Yp=cell(1,maxar); rss=zeros(1,maxar); BIC=zeros(1,maxar);

for order=minar:maxar
    BIC(order)=log(models{order}.noisevariance)+...
         (size(models{order}.parametervector,1)*log(N))/N;
end

for order=minar:maxar
    Yp{order}=predict(models{order},DATA,1);
    rss(order)=100*(norm(DATA.outputdata-Yp{order}.outputdata)^2)/(norm(DATA.outputdata)^2);
end

figure (4)
subplot(2,1,1),plot(minar:maxar,BIC(minar:maxar),'-o')
xlim([minar maxar])
title('Sensor 2')
ylabel('BIC')
xlabel('Model order')

subplot(2,1,2),plot(minar:maxar,rss(minar:maxar),'-o')
xlim([minar maxar])

ylabel('RSS/SSS (%)')
xlabel('Model order')

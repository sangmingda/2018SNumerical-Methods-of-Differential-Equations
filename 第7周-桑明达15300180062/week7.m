clear all;close all;clc

%% P118 1 $\frac{\mathrm{d}x}{\mathrm{d}t}=\lambda \left ( -u+cos\left ( t \right ) \right )$
clear all;close all;clc

u0=[0,1];
lambda=[1,1e1,1e2,1e3];
C=u0'-lambda.^2./(1+lambda.^2); %2*4矩阵
title0={'u(0)=0,\lambda =1' 'u(0)=0,\lambda =10' 'u(0)=0,\lambda =100' 'u(0)=0,\lambda =1000'
    'u(0)=1,\lambda =1' 'u(0)=1,\lambda =10' 'u(0)=1,\lambda =100' 'u(0)=1,\lambda =1000'};
    
 for i=1:length(u0)
 for j=1:length(lambda)

u_0=u0(1,i);
lambda0=lambda(1,j);
C0=C(i,j);
Ut=@(t) lambda0./(1+lambda0.^2).*(sin(t)+lambda0.*cos(t))+C0.*exp(-lambda0.*t);
F1=@(u,t,dt) u+dt.*lambda0.*(-u+cos(t));
F2=@(u,t,dt) (u+dt.*lambda0.*cos(t))./(1+dt.*lambda0);

a=0;b=10;n=101;%积分区间
dt=(b-a)./(n-1);
t=linspace(a,b,n);
y=Ut(t);
% plot(x,y)

y1=u_0;
y2=u_0;

for k=1:(n-1)
y1(k+1)=F1(y1(k),t(k),dt);%显式Euler格式
y2(k+1)=F2(y2(k),t(k+1),dt);%隐式Euler格式
end

figure(4.*i+j-4)
plot(t,y);hold on;
plot(t,y1,'^r');hold on;
plot(t,y2,'*b');
legend('精确值','显式Euler格式','隐式Euler格式');
title(title0{i,j});

 end
 end
 %%
 j=4;%\lambda=1000;
 for i=1:length(u0)

u_0=u0(1,i);
lambda0=lambda(1,j);
C0=C(i,j);

Ut=@(t) lambda0./(1+lambda0.^2).*(sin(t)+lambda0.*cos(t))+C0.*exp(-lambda0.*t);

F=@(u,t) lambda0.*(-u+cos(t));
Adams=@(u,t,dt) dt./12.*(23.*F(u(3),t(3))-16.*F(u(2),t(2))+5.*F(u(1),t(1)));%三步三阶Adams-Bashforth方法
Gear=@(u,t,dt) (-3.*u(3)+3/2.*u(2)-1/3.*u(1)-dt.*lambda0.*cos(t))./(-11/6-dt.*lambda0);%三阶Gear格式

a=0;b=10;n=101;%积分区间
dt=(b-a)./(n-1);
t=linspace(a,b,n);
y=Ut(t);
% plot(x,y)

y1=y(1:3);
y2=y(1:3);

for k=3:(n-1)
y1(k+1)=y1(k)+Adams(y1((k-2):k),t((k-2):k),dt);%三步三阶Adams-Bashforth方法
y2(k+1)=Gear(y2((k-2):k),t(k+1),dt);%三阶Gear格式
end

figure(2*i+7)
plot(t,y);hold on;
%plot(t,y1,'^r');hold on;
plot(t,y2,'*b');
legend('精确值','三阶Gear格式');
title(title0{2*j+i-2});

figure(2*i+8)
semilogy(t,y);hold on;
semilogy(t,y1,'^r');hold on;
semilogy(t,y2,'*b');
legend('精确值','三步三阶Adams-Bashforth方法','三阶Gear格式');
title(title0{2*j+i-2});

 
 end
 
 
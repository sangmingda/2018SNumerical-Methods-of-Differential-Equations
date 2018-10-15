%% week4.1 P69 1
clear all;close all;clc




%% week4.2 P69 2 Euler格式计算$\frac{du}{dt} = au$
clear all;close all;clc

T=1;%区间长度
deltat=0.04;%步长
u0=1;%初始值

% 求解$\frac{du}{dt} = f$
% f所求方程
% T区间长度
% deltat步长
% x0初始值
%y1,y2,y3,y4分别是显式，隐式，改进，修正的Euler格式结果

a=-2;%方程系数
f=@(t,x) a*x;

n=floor(T(1)/deltat(1));
y1=u0;
y2=u0;
y3=u0;
y4=u0;
t=u0+deltat.*(0:1:n);

for i=1:n
y1(i+1)=y1(i)+deltat(1).*f(t,y1(i));
y2(i+1)=y2(i)./(1-deltat(1).*a);%需要根据f修改
y3(i+1)=(1+deltat(1).*a./2).*y3(i)./(1-deltat(1).*a./2);%需要根据f修改
y4(i+1)=(1+deltat(1).*a.*(1+deltat(1).*a./2)).*y4(i);%需要根据f修改
end
y=u0.*exp(a.*(t-1));
plot(t,y);
hold on
scatter(t,y1,'ro');
hold on
scatter(t,y2,'g+');
hold on
scatter(t,y3,'b*');
hold on
scatter(t,y4,'m.');
legend('真实解','显式Euler格式','隐式Euler格式','改进Euler格式','修正Euler格式');
%% week4.3 P73 1 改进，修正的Euler格式的稳定性分析和绝对稳定区间
clear all;close all;clc



%% week4.4 P74 1 Taylor级数计算
clear all;close all;clc

t0=0;
t1=8;%区间长度
deltat=0.2;%步长

f=@(x) x-x.^2;
F=@(x) (x-x.^2).*(1-2.*x);
G=@(x) (x-x.^2).^2.*(-2);
fu=@(x) (1-2.*x);

n=floor((t1-t0)/deltat(1));
t=t0+deltat.*(0:1:n);

u0=0.5;%初始值
y1=u0;
y2=u0;
y3=u0;

for i=1:n
y1(i+1)=y1(i)+deltat.*f(y1(i));
y2(i+1)=y2(i)+deltat.*(f(y2(i))+0.5.*deltat.*F(y2(i)));
y3(i+1)=y3(i)+deltat.*(f(y3(i))+0.5.*deltat.*F(y3(i))+deltat.^2.*(G(y3(i))+fu(y3(i)).*F(y3(i)))./6);
end
y=1./((1./u0-1).*exp(-t)+1);

plot(t,y);
axis([0 8 0.5 1.5]); 
hold on
scatter(t,y1,'ro');
hold on
scatter(t,y2,'gd');
hold on
scatter(t,y3,'b*');
hold on

u0=1.5;%初始值
y4=u0;
y5=u0;
y6=u0;

for i=1:n
y4(i+1)=y4(i)+deltat.*f(y4(i));
y5(i+1)=y5(i)+deltat.*(f(y5(i))+0.5.*deltat.*F(y5(i)));
y6(i+1)=y6(i)+deltat.*(f(y6(i))+0.5.*deltat.*F(y6(i))+deltat.^2.*(G(y6(i))+fu(y6(i)).*F(y6(i)))./6);
end
y0=1./((1./u0-1).*exp(-t)+1);

plot(t,y0);
hold on
scatter(t,y4,'ro');
hold on
scatter(t,y5,'gd');
hold on
scatter(t,y6,'b*');

legend('真实解','显式Euler格式','q=2 ','q=3','真实解','显式Euler格式','q=2 ','q=3');



%% week4.5 P79 2 例2.3.1
clear all;close all;clc
for j=1:15
t0=1;
t1=2;%区间长度
deltat=2^(-j);%步长

f=@(x) -2.*x;

n=floor((t1-t0)/deltat(1));
t=t0+deltat.*(0:1:n);

u0=0.5;%初始值
y1=u0;
y2=u0;
y3=u0;
y4=u0;
for i=1:n
y1(i+1)=y1(i)+deltat.*f(y1(i));
y2(i+1)=y2(i)+deltat./2.*(f(y2(i))+f(y2(i)+deltat.*f(y2(i))));
y3(i+1)=y3(i)+deltat./4.*(f(y3(i))+3.*f(y3(i)+2.*deltat./3.*f(y3(i)+deltat./3.*f(y3(i)))));
y4(i+1)=y4(i)+deltat./6.*(f(y4(i))+2.*f(y4(i)+deltat./2.*f(y4(i)))+2.*f(y4(i)+deltat./2.*f(y4(i)+deltat./2.*f(y4(i))))+f(y4(i)+deltat.*f(y4(i)+deltat./2.*f(y4(i)+deltat./2.*f(y4(i))))));
end
y=u0.*exp(-2.*(t-t0));
R1(j)=abs(y1(end)-y(end));
R2(j)=abs(y2(end)-y(end));
R3(j)=abs(y3(end)-y(end));
R4(j)=abs(y4(end)-y(end));
end
figure(1)
plot(1:j,log(R1),'.');
hold on
plot(1:j,log(R2),'ro');
hold on
plot(1:j,log(R3),'gd');
hold on
plot(1:j,log(R4),'b*');
legend('显式Euler方法','二阶方法 ','三阶方法','四阶方法');

figure(2)
plot(t,y);
hold on
scatter(t,y1,'ro');
hold on
scatter(t,y2,'gd');
hold on
scatter(t,y3,'b*');
hold on
scatter(t,y4,'m.');
legend('真实解','显式Euler方法','二阶方法 ','三阶方法','四阶方法');

%% week4.6 P79 3 例2.2.2
clear all;close all;clc

t0=1;
t1=8;
u0=0.5;%初始值

f=@(x) x-x.^2;
F=@(x) 1./((1./u0-1).*exp(-x+t0)+1);

for j=1:15

deltat=2^(-j);%步长

n=floor((t1-t0)/deltat(1));
t=t0+deltat.*(0:1:n);

y1=u0;
y2=u0;
y3=u0;
y4=u0;
for i=1:n
y1(i+1)=y1(i)+deltat.*f(y1(i));
y2(i+1)=y2(i)+deltat./2.*(f(y2(i))+f(y2(i)+deltat.*f(y2(i))));
y3(i+1)=y3(i)+deltat./4.*(f(y3(i))+3.*f(y3(i)+2.*deltat./3.*f(y3(i)+deltat./3.*f(y3(i)))));
y4(i+1)=y4(i)+deltat./6.*(f(y4(i))+2.*f(y4(i)+deltat./2.*f(y4(i)))+2.*f(y4(i)+deltat./2.*f(y4(i)+deltat./2.*f(y4(i))))+f(y4(i)+deltat.*f(y4(i)+deltat./2.*f(y4(i)+deltat./2.*f(y4(i))))));
end
y=F(t);
R1(j)=abs(y1(end)-y(end));
R2(j)=abs(y2(end)-y(end));
R3(j)=abs(y3(end)-y(end));
R4(j)=abs(y4(end)-y(end));
end
figure(1)
plot(1:j,log(R1),'.');
hold on
plot(1:j,log(R2),'ro');
hold on
plot(1:j,log(R3),'gd');
hold on
plot(1:j,log(R4),'b*');
legend('显式Euler方法','二阶方法 ','三阶方法','四阶方法');

figure(2)
plot(t,y);
hold on
scatter(t,y1,'ro');
hold on
scatter(t,y2,'gd');
hold on
scatter(t,y3,'b*');
hold on
scatter(t,y4,'m.');
legend('真实解','显式Euler方法','二阶方法 ','三阶方法','四阶方法');
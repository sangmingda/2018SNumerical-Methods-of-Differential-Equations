clear all;close all;clc

%% P195 1 抛物型方程
clear all;close all;clc

a=1;

xnum=100;
%x=linspace(0,1,xnum+1);
u0=zeros(xnum+1,1);
f=@(x,t) sin(pi*x)+sin(2*pi*x);

t=5;
tnum=100;
theta=0;

%u=theta_func(u0,a,f,t,tnum,theta);
u=Richardson_func(u0,a,f,t,tnum);

subplot(1,3,1)
contour(u,'ShowText','on')
set(gca,'YTick',0:0.25*tnum:tnum);
set(gca,'XTick',0:0.25*xnum:xnum);
set(gca,'YTickLabel',{'0','1/4','1/2','3/4','1'});
set(gca,'XTickLabel',{'0',num2str(1/4*t),num2str(1/2*t),num2str(3/4*t),num2str(1*t)});
axis([0 tnum 0 xnum])
xlabel('时间t')
ylabel('距离x')
%title(['\theta =',num2str(theta)])
title('Richardson格式')
%误差
U0=@(x,t) (1-exp(-pi^2.*t))./(pi^2).*sin(pi*x)+(1-exp(-4*pi^2.*t))./(4*pi^2).*sin(2*pi*x);
for i=1:tnum+1
    for j=1:xnum+1
U_0(j,i)=U0((j-1)/(xnum),(i-1)*t*(tnum));
    end
end

subplot(1,3,2)
contour(U_0,'ShowText','on')
set(gca,'YTick',0:0.25*tnum:tnum);
set(gca,'XTick',0:0.25*xnum:xnum);
set(gca,'YTickLabel',{'0','1/4','1/2','3/4','1'});
set(gca,'XTickLabel',{'0',num2str(1/4*t),num2str(1/2*t),num2str(3/4*t),num2str(1*t)});
axis([0 tnum 0 xnum])
xlabel('时间t')
ylabel('距离x')
title('真实解')

error=abs(u-U_0)./abs(U_0);
error=log(error);

subplot(1,3,3)
contour(error,'ShowText','on')
set(gca,'YTick',0:0.25*tnum:tnum);
set(gca,'XTick',0:0.25*xnum:xnum);
set(gca,'YTickLabel',{'0','1/4','1/2','3/4','1'});
set(gca,'XTickLabel',{'0',num2str(1/4*t),num2str(1/2*t),num2str(3/4*t),num2str(1*t)});
axis([0 tnum 0 xnum])
xlabel('时间t')
ylabel('距离x')
title('误差(取对数)')


%% P195 1 抛物型方程
clear all;close all;clc

a=1;

xnum=100;
x=linspace(0,1,xnum+1);
u0=[0;rand(xnum-1,1);0];
f=@(x,t) sin(pi*x.*exp(-t.^2));

t=5;
tnum=100;

u1=theta_func(u0,a,f,t,tnum,0);
u2=theta_func(u0,a,f,t,tnum,0.5);
u3=theta_func(u0,a,f,t,tnum,1);
u4=Richardson_func(u0,a,f,t,tnum);

subplot(1,4,1)
contour(u1,'ShowText','on')
set(gca,'YTick',0:0.25*tnum:tnum);
set(gca,'XTick',0:0.25*xnum:xnum);
set(gca,'YTickLabel',{'0','1/4','1/2','3/4','1'});
set(gca,'XTickLabel',{'0',num2str(1/4*t),num2str(1/2*t),num2str(3/4*t),num2str(1*t)});
axis([0 tnum 0 xnum])
xlabel('时间t')
ylabel('距离x')
title('\theta = 0')

subplot(1,4,2)
contour(u3,'ShowText','on')
set(gca,'YTick',0:0.25*tnum:tnum);
set(gca,'XTick',0:0.25*xnum:xnum);
set(gca,'YTickLabel',{'0','1/4','1/2','3/4','1'});
set(gca,'XTickLabel',{'0',num2str(1/4*t),num2str(1/2*t),num2str(3/4*t),num2str(1*t)});
axis([0 tnum 0 xnum])
xlabel('时间t')
ylabel('距离x')
title('\theta = 0.5')

subplot(1,4,3)
contour(u3,'ShowText','on')
set(gca,'YTick',0:0.25*tnum:tnum);
set(gca,'XTick',0:0.25*xnum:xnum);
set(gca,'YTickLabel',{'0','1/4','1/2','3/4','1'});
set(gca,'XTickLabel',{'0',num2str(1/4*t),num2str(1/2*t),num2str(3/4*t),num2str(1*t)});
axis([0 tnum 0 xnum])
xlabel('时间t')
ylabel('距离x')
title('\theta = 1')

subplot(1,4,4)
contour(u4,'ShowText','on')
set(gca,'YTick',0:0.25*tnum:tnum);
set(gca,'XTick',0:0.25*xnum:xnum);
set(gca,'YTickLabel',{'0','1/4','1/2','3/4','1'});
set(gca,'XTickLabel',{'0',num2str(1/4*t),num2str(1/2*t),num2str(3/4*t),num2str(1*t)});
axis([0 tnum 0 xnum])
xlabel('时间t')
ylabel('距离x')
title('Richardson格式')



%%
function u=theta_func(u0,a,f,t,tnum,theta)

u0=reshape(u0,[],1);
xnum=size(u0,1)-1;
h=1/(xnum);
tau=abs(t)./(tnum);
Delta_h=(diag(ones(1,xnum-2),1)-2.*diag(ones(1,xnum-1),0)+diag(ones(1,xnum-2),-1))./h^2;
x0=linspace(0,1,xnum+1);
t0=linspace(0,t,tnum+1);
u=[u0 zeros(xnum+1,tnum)];
for i=2:tnum+1
    A=diag(ones(1,xnum-1),0)-a.*tau.*(theta).*Delta_h;
     B=diag(ones(1,xnum-1),0)+a.*tau.*(1-theta).*Delta_h;
     C=tau.*(theta.*f(x0(2:xnum),t0(i))+(1-theta).*f(x0(2:xnum),t0(i)));
    u(2:xnum,i)=A\(B*u(2:xnum,i-1)+C');  
end
end

function u=theta0_func(u0,a,f,t,tnum)

u0=reshape(u0,[],1);
xnum=size(u0,1)-1;
h=1/(xnum);
tau=abs(t)./(tnum);
Delta_h=(diag(ones(1,xnum-2),1)-2.*diag(ones(1,xnum-1),0)+diag(ones(1,xnum-2),-1))./h^2;
x0=linspace(0,1,xnum+1);
t0=linspace(0,t,tnum+1);
u=[u0 zeros(xnum+1,tnum)];
for i=2:tnum+1
   
     B=diag(ones(1,xnum-1),0)+a.*tau.*Delta_h;
     C=tau.*(f(x0(2:xnum),t0(i)));
    u(2:xnum,i)=(B*u(2:xnum,i-1)+C');  
end
end



function u=Richardson_func(u0,a,f,t,tnum)

u0=reshape(u0,[],1);
xnum=size(u0,1)-1;
h=1/(xnum);
tau=abs(t)./(tnum);
Delta_h=(diag(ones(1,xnum-2),1)-2.*diag(ones(1,xnum-1),0)+diag(ones(1,xnum-2),-1))./h^2;
x0=linspace(0,1,xnum+1);
t0=linspace(0,t,tnum+1);
u=[u0 zeros(xnum+1,tnum)];
u(2:xnum,2)=2*tau.*(f(x0(2:xnum),t0(1)))'+2*a.*tau.*Delta_h*u(2:xnum,1);
for i=3:tnum+1   
     B=2*a.*tau.*Delta_h;
     C=2*tau.*(f(x0(2:xnum),t0(i)));
    u(2:xnum,i)=(u(2:xnum,i-2)+B*u(2:xnum,i-1)+C');  
end
end
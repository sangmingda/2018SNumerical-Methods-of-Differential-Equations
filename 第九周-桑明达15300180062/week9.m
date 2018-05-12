clear all;close all;clc

%% P143 3 三点差分及精度
clear all;close all;clc

a=1;
b=1;
c=1;
f=@(x) 1.*x.^0;
u0=0;
%uN=0:
x0=0;xN=1;%求解区间

n=100;
h=(xN-x0)./n;

A=(-2.*diag(ones(1,n-1),0)+1.*diag(ones(1,n-2),1)+1.*diag(ones(1,n-2),-1));
A(n-1,n-1)=A(n-1,n-1)+1./(h+1);
A=-1./(h.^2).*A;

B=(-1.*diag(ones(1,n-2),1)+1.*diag(ones(1,n-2),-1));
B(n-1,n-1)=B(n-1,n-1)+1./(h+1);
B=b./(2.*a.*h).*B;

C=c./a.*eye(n-1,n-1);

x=linspace(x0,xN,n+1);
fx=f(x);
F=fx(1,2:(n))';
F(1)=F(1)+u0./(h.^2)+b.*u0./(2.*a.*h);


u=(A+B+C)\F;

uN=u(end)./(h+1);
U=[u0,u',uN];
plot(x,U)








clear all;close all;clc

m=10;n=10;
a=0;b=1;

for i=1:m
    for j=1:n
        f=@(x) x.^i;
        I(i,j)=newton_cotes(f,a,b,j-1);        
    end
end

function I=newton_cotes(f,a,b,n)
% f是被积函数；
% a,b是区间范围；
% n是Newton公式阶数；
n=floor(n);
if n<2
    I=f(0.5*a+0.5*b);
else
    x=linspace(a,b,n+1);
    y=f(x);
    y=y';

c1=1:n;
c2=1:n;
c3=c1.^(c2');
c4=0;
for i=1:n
    c4(i,1)=n.^(i)./(i+1);
end
%c_nk=(c3)^(-1)*c4;
c_nk=(c3)\c4;
c_nk=[c_nk(n);c_nk];
I=(b-a).*((c_nk')*y);
end
end
clear all;close all;clc

m=10;n=10;
a=0;b=1;

for i=1:m
    for j=1:n
        f=@(x) x.^i;
        I(i,j)=gauss(f,a,b,j-1);        
    end
end

function G=gauss(f,a,b,n)
% f是被积函数；
% a,b是区间范围；
% n是Guass公式阶数；
n=floor(n);
if n<1
    G=f(0.5*a+0.5*b);
else
    %生成求积点Gaussx和权w
    %n=1;
    B=.5./sqrt(1-(2*(1:n)).^(-2));
    [V,Lambda]=eig(diag(B,1)+diag(B,-1));
    [x,Isort]=sort(diag(Lambda));
    w=2*V(1,Isort).^2;
    x=(b-a)./2.*x+(b+a)./2;

    %计算积分值
    G=0;
for k=1:length(x)
    G=G+(b-a)./2.*w(k).*f(x(k));
end
end
end
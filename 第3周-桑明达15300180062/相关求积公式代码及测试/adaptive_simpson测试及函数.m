clear all;close all;clc

f=@(x) (x)^0.5;
a=0;
b=1;
tol=15*eps;
bottom=10;

I=adaptive_simpson(f,a,b,tol,bottom)


function I=adaptive_simpson(f,a,b,tol,bottom)
% f是被积函数；
% （a，b）是求积区间；
% tol是误差限制；
% bottom是最大计算次数。
a=a(1);b=b(1);
f0=f(a);f1=f(0.75*a+0.25*b);f2=f(0.5*a+0.5*b);f3=f(0.25*a+0.75*b);f4=f(b);
S0=(b-a)*(f0+4*f2+f4)/6;
S1=(b-a)*(f0+4*f1+2*f2+4*f3+f4)/12;

if abs(S0-S1)<tol 
    I=S1+(S1-S0)/15;
else    
    I=adaptive_simpson(f,a,0.5*a+0.5*b,tol,bottom-1)+adaptive_simpson(f,0.5*a+0.5*b,b,tol,bottom-1);    
end

end
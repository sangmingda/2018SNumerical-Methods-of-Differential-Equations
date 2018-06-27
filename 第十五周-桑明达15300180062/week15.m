clear all;close all;clc

%% P254 3 用线性有限元求解混合边值问题
clear all;close all;clc

for m=1:6

h=2^(-m);
n=2^m;
x=linspace(0,1,n+1)';

u0=1;%第一类边界条件
u_1=-pi/2;%第二类边界条件

l_0=2*(1./h-h./3);
l_1=-1./h-h./6;
L=diag(l_0*ones(1,n-1),0)+diag(l_1*ones(1,n-2),1)+diag(l_1*ones(1,n-2),-1);
L(n-1,n-1)=L(n-1,n-1)+l_1;

for i=2:n   %f1 to fN-1
f1=@(xx) (pi^2/4+1)*cos(pi*xx/2).*(xx-x(i-1))/h;
f2=@(xx) (pi^2/4+1)*cos(pi*xx/2).*(x(i+1)-xx)/h;
F(i-1,1) =integral(f1,x(i-1),x(i))+integral(f2,x(i),x(i+1));
end
F(1)=F(1)-u0*l_1;
F(n-1)=F(n-1)-h*u_1*l_0;

uu=L\F;
uu=[u0;uu;h*u_1+uu(end)];
% uux=(uu(1:n-1,1)-[u0;uu(1:n-2,1)])./h;

U=@(x) cos(pi*x/2);
% Ux=@(x) -pi/2*sin(pi*x/2);
y=U(x);
%yx=Ux(x);

figure(m)
plot(x,uu,'*',x,y,'.-')
thetitle=['h=2^',num2str(m)];
title(thetitle)
legend('计算解','精确解')
thefilename=['week15_1_',num2str(m),'.eps'];
saveas(gcf,thefilename)

fprintf('\\begin{figure}[H]\n')
fprintf('\\includegraphics[width=1\\linewidth]{')
fprintf(thefilename)
figname=['}\n\\label{Fig:',num2str(m),'}\n\\end{figure}\n\n'];
fprintf(figname)
end
% 	\begin{figure}[H]
% 		\includegraphics[width=1\linewidth]{week15_1_1.eps}
% 		\label{Fig:1}
% 	\end{figure}

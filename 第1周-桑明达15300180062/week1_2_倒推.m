%%
clear all;close all;clc

for j=1:100
    n=100*j;
u(n,1)=0.1;
u(n-1,1)=0.1;
for i=(n-2):-1:1
    u(i,1)=(3*u(i+1,1)-u(i+2,1))./2;
end
tol(j,1)=u(1)-0.1;
end
fprintf('%.40s\n',tol);
fprintf('%.40s\n',u(1));
fprintf('%.40s\n',2*tol(1));
%%
clear all;close all;clc
A=[1.5 -0.5
    1 0];
[V,D]=eig(A)
V^-1
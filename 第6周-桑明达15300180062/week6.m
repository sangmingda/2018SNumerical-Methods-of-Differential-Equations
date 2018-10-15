clear all;close all;clc

%% P93 2 不同的初始选取对精度的影响
clear all;close all;clc

a = -2;
f = @(t,u) a*u;
t0 = 0;
u0 = 1;
T =0.31;
dt = 1e-1;

u_exm=exp(-2.*(0:0.1:1));

ilist = [1,3,5,7];
method_list = {'Euler-Explicit', ...
               'Euler-Modified', 'Runge-Kutta2', 'Heun2', ...
               'Kutta3', 'Heun3', ...
               'Runge-Kutta4', 'Kutta4'};
u1=zeros(4,4);
for i=1:4
    method = method_list{ilist(i)};
    
    u = Ch2_RungeKutta(f,t0,u0,dt,T,method);
    u1(i,:)=u1(i,:)+u';
    
end
u2=exp(-2.*(0:0.1:0.3));
u0=[u2
    u1];


ti = 0;
T =1.01;
k=4;
u1=zeros(5,11);
u2=zeros(5,11);
for i=1:5
    
    ui = u0(i,:);
    u = Ch2_AdamsBashforth(f,ti,ui,dt,T,k);
    u1(i,:)=u1(i,:)+u';
    u = Ch2_Geer(f,ti,ui,dt,T,k,a);
    u2(i,:)=u2(i,:)+u'; 
end
err1=abs(u1-u_exm);
err2=abs(u2-u_exm);
x=0:0.1:1;
markers = {'r+-','gd-','b*-','co-','y--o'};
figure(1)
for i=1:5
semilogy(x,err1(i,:),markers{i}); hold on;

end
names = {'精确的初始值','EulerExplicit','Runge-Kutta2','Kutta3','Runge-Kutta4'};
legend(names,'Location','SouthEast');
title('四阶Adams格式精度测试（误差曲线）');

figure(2)
for i=1:5
semilogy(x,err2(i,:),markers{i}); hold on;

end
names = {'精确的初始值','EulerExplicit','Runge-Kutta2','Kutta3','Runge-Kutta4'};
legend(names,'Location','SouthEast');
title('四阶Geer格式精度测试（误差曲线）');
%% 函数
function u = Ch2_RungeKutta(f,ti,ui,dt,T,method)
%--------------------------------------------------------------------------
% u = Ch2_RungeKutta(f,ti,ui,dt,T,method)
% 
% Runge-Kutta method for solving ODE
%  method: 'Euler-Explicit';
%          'Midpoint'='Euler-Modified', 'Runge-Kutta2', 'Heun2';
%          'Kutta3', 'Heun3';
%          'Runge-Kutta4', 'Kutta4';
% 
% Example by Wenbin Chen & Xinming Wu
%               2012-10-01
%--------------------------------------------------------------------------

N = floor((T-ti)/dt);
u = zeros(N+1,1);
u(1) = ui;

for n=1:N
    t0 = ti + dt*(n-1);
    u0 = u(n);
    
    switch method
        case 'Euler-Explicit'
            k1 = f(t0,u0);
            u1 = u0 + dt*k1;
            
        case {'Midpoint', 'Euler-Modified'}
            k1 = f(t0,u0);
            k2 = f(t0+dt/2,u0+dt/2*k1);
            u1 = u0 + dt*k2;
            
        case 'Runge-Kutta2'
            k1 = f(t0,u0);
            k2 = f(t0+dt,u0+dt*k1);
            u1 = u0 + dt*(k1+k2)/2;
            
        case 'Heun2'
            k1 = f(t0,u0);
            k2 = f(t0+2*dt/3,u0+2*dt/3*k1);
            u1 = u0 + dt*(k1+3*k2)/4;
            
        case 'Kutta3'
            k1 = f(t0,u0);
            k2 = f(t0+dt/2,u0+dt/2*k1);
            k3 = f(t0+dt,u0-dt*k1+2*dt*k2);
            u1 = u0 + dt*(k1+4*k2+k3)/6;
            
        case 'Heun3'
            k1 = f(t0,u0);
            k2 = f(t0+dt/3,u0+dt/3*k1);
            k3 = f(t0+2*dt/3,u0+2*dt/3*k2);
            u1 = u0 + dt*(k1+3*k3)/4;
            
        case 'Runge-Kutta4'
            k1 = f(t0,u0);
            k2 = f(t0+dt/2,u0+dt/2*k1);
            k3 = f(t0+dt/2,u0+dt/2*k2);
            k4 = f(t0+dt,u0+dt*k3);
            u1 = u0 + dt*(k1+2*k2+2*k3+k4)/6;
            
        case 'Kutta4'
            k1 = f(t0,u0);
            k2 = f(t0+dt/3,u0+dt/3*k1);
            k3 = f(t0+2*dt/3,u0-dt/3*k1+dt*k2);
            k4 = f(t0+dt,u0+dt*k1-dt*k2+dt*k3);
            u1 = u0 + dt*(k1+3*k2+3*k3+k4)/8;
    end
    
    u(n+1) = u1;
end
end

function u = Ch2_AdamsBashforth(f,ti,ui,dt,T,k)
%--------------------------------------------------------------------------
% u = Ch2_AdamsBashforth(f,ti,ui,dt,T,k)
% 
% Adams-Bashforth method for solving ODE
%  k: order of the method
% 
% Example by Wenbin Chen & Xinming Wu
%               2012-10-01
%--------------------------------------------------------------------------

N = floor((T-ti)/dt);
u = zeros(N+1,1);
fval = zeros(N+1,1);

t = (ti:dt:T).';
u(1:k) = ui;
fval(1:k) = f(t(1:k),u(1:k));

for n=k:N
    u0 = u(n);
    t1 = t(n+1);
    
    switch k
        case 1
            u1 = u0 + dt*fval(n);
            
        case 2
            u1 = u0 + dt*[3,-1]/2*fval(n:-1:n-1);
            
        case 3
            u1 = u0 + dt*[23,-16,5]/12*fval(n:-1:n-2);
            
        case 4
            u1 = u0 + dt*[55,-59,37,-9]/24*fval(n:-1:n-3);
            
    end
    
    u(n+1) = u1;
    fval(n+1) = f(t1,u1);
end
end

function u = Ch2_Geer(f,ti,ui,dt,T,k,a)

N = floor((T-ti)/dt);
u = zeros(N+1,1);
fval = zeros(N+1,1);

t = (ti:dt:T).';
u(1:k) = ui;
fval(1:k) = f(t(1:k),u(1:k));

for n=k:N
    u0 = u(n);
    t1 = t(n+1);
    
    switch k
        case 1
            u1 =(u(n))./(-a*dt+1);
            
        case 2
            u1 =(+[2,-1./2]*u((n):-1:(n-1)))./(-a*dt+3/2);
            
        case 3
            u1 =(+[3,-3/2,1/3]*u((n):-1:(n-2)))./(-a*dt+11/6);
            
        case 4
            u1 =(-[-4,3,-4/3,1/4]*u((n):-1:(n-3)))./(-a*dt+25/12);
            
        case 5
            u1 =(-[-5,5,-10/3,5/4,-1/5]*u((n):-1:(n-4)))./(-a*dt+137/60);
            
        case 6
            u1 =(-[-6,15/2,-20/3,15/4,-6/5,1/6]*u((n):-1:(n-5)))./(-a*dt+147/60);
            
    end
    
    u(n+1) = u1;
    fval(n+1) = f(t1,u1);
end
end
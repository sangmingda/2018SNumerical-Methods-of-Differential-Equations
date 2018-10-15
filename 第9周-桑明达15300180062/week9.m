clear all;close all;clc

%% P143 3 三点差分及精度
clear all;close all;clc

B0=[0,0,10,-10];
C0=[0,10,0,0];
markers = {'r+-','gd-','b*-','co-','y--o'};
names={'b=0,c=0','b=0,c=10','b=10,c=0','b=-10,c=0'};
for i=1:4


a=1;
b=B0(i);
c=C0(i);
f=@(x) 1.*x.^0;
u0=0;
%uN=0:
x0=0;xN=1;%求解区间

n=100;
h=(xN-x0)./n;

A=(-2.*diag(ones(1,n-1),0)+1.*diag(ones(1,n-2),1)+1.*diag(ones(1,n-2),-1));
A(n-1,n-1)=A(n-1,n-1)+1./(h+1);
A=-a./(h.^2).*A;

B=(1.*diag(ones(1,n-2),1)-1.*diag(ones(1,n-2),-1));
B(n-1,n-1)=B(n-1,n-1)+1./(h+1);
B=b./(2.*h).*B;

C=c.*eye(n-1,n-1);

x=linspace(x0,xN,n+1);
fx=f(x);
F=fx(1,2:(n))';
F(1)=F(1)+u0./(h.^2)+b.*u0./(2.*a.*h);


u=(A+B+C)\F;

uN=u(end)./(h+1);
U(i,:)=[u0,u',uN];%三点差分结果
figure(1)
hold on
plot(x,U(i,:),markers{i})

lambda1=(-b-sqrt(b.^2+4.*a.*c))./(-2.*a);
lambda2=(-b+sqrt(b.^2+4.*a.*c))./(-2.*a);
alpha1=-(exp(lambda2).*(1+lambda2)-1)./(c.*(exp(lambda2).*(1+lambda2)-exp(lambda1).*(1+lambda1)));
alpha2=-(1-exp(lambda1).*(1+lambda1))./(c.*(exp(lambda2).*(1+lambda2)-exp(lambda1).*(1+lambda1)));
if c==0
        if b==0
            Ux=@(x) x.*(3-2.*x)./(4.*a);
        else
            Ux=@(x) -2.*(exp(b.*x./a)-1)./(b.*((b./a+1).*exp(b./a)-1))+x./b;
        end
else
    Ux=@(x) alpha1.*exp(lambda1.*x)+alpha2.*exp(lambda2.*x)+1./c;
end

U0(i,:)=Ux(x);

figure(2)
hold on
loglog(x,U0(i,:),markers{i})

error(i)=sum(abs(U0(i,:)-U(i,:)))./sum(abs(U0(i,:)));
fprintf(names{i})
fprintf('时误差率是%s\n',error(i))
end
figure(1)
title('三点差分格式离散求解');
legend(names);


figure(2)
title('精确解');
%names={'b=0,c=0','b=0,c=10','b=10,c=0','b=-10,c=0'};
legend(names);





%% P146 图3.5 标准三点差分格式和四阶HOC格式
clear all;close all;clc

% Sec 3.1.4// Fig 3.5
a = 1; b = 1; c = 1;
ufun = @(x) sin(pi*x);
f = @(x) (a*pi^2+c)*sin(pi*x) + b*pi*cos(pi*x);

errs = zeros(11,2);
for i=1:11
    N = 2^(i+1);
    x = linspace(0,1,N+1).';
    uex = ufun(x);
    u = Ch3_FD1d(a,b,c,f,N);
    uhoc = Ch3_FD1dHOC(b,f,N);
    errs(i,1) = norm(u-uex)/norm(uex);
    errs(i,2) = norm(uhoc-uex)/norm(uex);
end
hs = 2.^(-(1:11)-1);
figure(1); clf; loglog(hs,errs(:,1),'r.-',hs,errs(:,2),'bo-');
legend('3-point-fd','HOC4','Location','NorthWest');
title('High Order Compact scheme');


function u = Ch3_FD1d(a,b,c,f,N,method)
%--------------------------------------------------------------------------
% u = Ch3_FD1d(a,b,c,f,N,method)
% 
% compute the numerical solution of the 1D example (h=1/N)
%   -a*u'' + b*u' + c*u = f  in (0,1)
%   u(0) = u(1) = 0
% 
% method for b: 'left'='backward', 'right'='forward', 'central'
% 
% Example by Wenbin Chen & Xinming Wu
%               2012-10-01
%--------------------------------------------------------------------------

if nargin<6, method = 'central'; end

h = 1/N;
e = ones(N-1,1);
A = spdiags([-e, 2*e, -e],-1:1,N-1,N-1)/h^2;
C = speye(N-1);

switch method
    case {'left','backward'} % b>0
        B = spdiags([-e, e],[-1,0],N-1,N-1)/h;
    case {'right','forward'} % b<0
        B = spdiags([-e, e],[0,1],N-1,N-1)/h;
    case {'center','central'}
        B = spdiags([-e, e],[-1,1],N-1,N-1)/(2*h);
end

if isa(f,'function_handle')
    x = (0:h:1)';
    F = f(x(2:end-1));
else
    F = f*ones(N-1,1);
end

u = (a*A + b*B + c*C) \ F;
u = [0; u; 0];
end

function u = Ch3_FD1dHOC(b,f,N)
%--------------------------------------------------------------------------
% u = Ch3_FD1dHOC(b,f,N)
% 
% compute the numerical solution of the 1D example (h=1/N)
%   -u'' + b*u' + u = f  in (0,1)
%   u(0) = u(1) = 0
% 
% Example by Wenbin Chen & Xinming Wu
%               2012-10-01
%--------------------------------------------------------------------------


h = 1/N;
e = ones(N-1,1);
A = spdiags([-e, 2*e, -e],-1:1,N-1,N-1)/h^2;
B = spdiags([-e, e],[-1,1],N-1,N-1)/(2*h); % 'central'
C = speye(N-1);

if isa(f,'function_handle')
    x = (0:h:1)';
    F = f(x(2:end-1));
    F0 = f(x(1)); F1 = f(x(end));
else
    F = f*ones(N-1,1);
    F0 = f; F1 = f;
end

h12 = h^2/12;
F = (C - b*h12*B - h12*A)*F;

% (-b*h^2/12)(-1/2h F0) + (-h^2/12)(-1/h^2 F0)
% (-b*h^2/12)(+1/2h F1) + (-h^2/12)(-1/h^2 F1)
F(1) = F(1) + b*h/24*F0 + 1/12*F0; 
F(end) = F(end) - b*h/24*F1 + 1/12*F1;

u = ((1+(b^2-1)*h12)*A + b*(1-h12)*B + C) \ F;
u = [0; u; 0];
end






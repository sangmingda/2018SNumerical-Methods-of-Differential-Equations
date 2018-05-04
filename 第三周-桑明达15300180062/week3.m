%% week3.1 P53 1
clear all;close all;clc
% 多种求积公式计算x^n 
% 求精确计算几次多项式

%% 1.Newton公式
clear all;close all;clc

m=10;%x^n的阶数
n=10;%Newton公式的阶数
a=0;b=1;
%计算积分
for i=1:m
    f=@(x) x.^i;
    I0(i,1)=integral(f,a,b);
end
for j=1:n
    for i=1:m        
        f=@(x) x.^i;
        I(i,j)=newton_cotes(f,a,b,j-1);
    end
end

%与精确值比较,以便于后面加粗
I_compare=zeros(m,n);
for j=1:n
    for i=1:m        
        if abs(I0(i,1)-I(i,j))<1e-10
        I_compare(i,j)=1;
        end
    end
end
%输出latex表格代码第一行
fprintf('\\begin{tabular}{')
for j=1:(n+2)
        fprintf('|c')
        if j==n+2
        fprintf('|')
        end
    end
fprintf('}\n')
%输出latex表格
fprintf('\n\\hline\n')
for i=1:m
    tex=['$x^{',num2str(i),'}$&'];
    fprintf(tex);
    fprintf('%.6f&',I0(i,1))
    for j=1:n
        %准确值加粗\textbf{}
        if I_compare(i,j)==1
            fprintf('{\\color{red}{\\textbf{')
        end
        %输出积分结果
        fprintf('%.6f',I(i,j))
        %准确值加粗\textbf{}
        %准确值变红{\color{red}{变红}}
        if I_compare(i,j)==1
            fprintf('}}}')
        end
        if j<n
        fprintf('&')
        end
    end
    %输出latex换行
    fprintf('\\\\\n\\hline\n')
end

%% Gauss求积公式
clear all;close all;clc

m=20;%x^n的阶数
n=10;%Newton公式的阶数
a=0;b=1;
%计算积分
for i=1:m
    f=@(x) x.^i;
    I0(i,1)=integral(f,a,b);
end
for i=1:m
    for j=1:n
        f=@(x) x.^i;
        I(i,j)=gauss(f,a,b,j-1);        
    end
end

%与精确值比较,以便于后面加粗
I_compare=zeros(m,n);
for j=1:n
    for i=1:m        
        if abs(I0(i,1)-I(i,j))<1e-10
        I_compare(i,j)=1;
        end
    end
end
%输出latex表格代码第一行
fprintf('\\begin{tabular}{')
for j=1:(n+2)
        fprintf('|c')
        if j==n+2
        fprintf('|')
        end
    end
fprintf('}\n')
%输出latex表格
fprintf('\n\\hline\n')
fprintf('$x^n$&精确值')
for i=1:n
    tex=['&',num2str(i),'点'];
        fprintf(tex);
        end
fprintf('\\\\\n\\hline\n')
for i=1:m
    tex=['$x^{',num2str(i),'}$&'];
    fprintf(tex);
    fprintf('%.6f&',I0(i,1))
    for j=1:n
        %准确值加粗\textbf{}
        if I_compare(i,j)==1
            fprintf('{\\color{red}{\\textbf{')
        end
        %输出积分结果
        fprintf('%.6f',I(i,j))
        %准确值加粗\textbf{}
        %准确值变红{\color{red}{变红}}
        if I_compare(i,j)==1
            fprintf('}}}')
        end
        if j<n
        fprintf('&')
        end
    end
    %输出latex换行
    fprintf('\\\\\n\\hline\n')
end


%% 用于调用的函数
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






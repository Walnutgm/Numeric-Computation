function [x,i,y]=PCG(A,b,epsilon)
%预处理共轭梯度法求解线性方程组Ax=b
%输入：系数矩阵A，向量b,收敛条件epsilon
%输出：数值解x，迭代步数i，残差向量的二范数y
[n,m]=size(A);
M=zeros(n,n); %M为jacobi预处理矩阵
for i=1:n
    M(i,i)=1/A(i,i);
end
x0=zeros(m,1);
r0=b-A*x0;
z0=M*r0;
u0=z0;
max=100;
%y=zeros(max,1);
for i=1:max
    alpha=(r0'*z0)/(u0'*A*u0);
    x=x0+alpha*u0;
    r=r0-alpha*A*u0;
    %y(i)=norm(r);
    z=M*r;
    beta=(z'*r)/(z0'*r0);
    u=z+beta*u0;
    %判断收敛条件1
%     if norm(x-x0)<=10^(-8)
%         break
%     end
    %判断收敛条件2
    if norm(r)<=epsilon
        break
    end
    x0=x;
    r0=r;
    z0=z;
    u0=u;
end
y=norm(b-A*x);
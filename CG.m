function [x,i,y]=CG(A,b,epsilon)
%共轭梯度法求解线性方程组Ax=b
%输入：系数矩阵A，向量b,收敛条件epsilon
%输出：数值解x，迭代步数i，残差向量的二范数y
[n,m]=size(A);
x0=zeros(m,1);
r0=b-A*x0;
u0=r0;
max=1000;
%y=zeros(max,1);
for i=1:max
    alpha=(r0'*r0)/(u0'*A*u0);
    x=x0+alpha*u0;
    r=r0-alpha*A*u0;
    %y(i)=norm(r);
    beta=(r'*r)/(r0'*r0);
    u=r+beta*u0;
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
    u0=u;
end
y=norm(b-A*x);

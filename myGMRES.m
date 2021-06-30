function [x,r]=myGMRES(A,b)
%使用gmres法求解线性方程组
%输入：系数矩阵A，右端向量b
%输出：求得的数值解x，所得解的残差的二范数r
[m, ~] = size(A);
x0=zeros(m,1);
H = zeros(m+1,m);
V = zeros(m,m+1);   %A*V=V*H
r0 = b-A*x0;
beta=norm(r0);
V(:,1) = r0./beta;
for j = 1:m
    w = A*V(:,j);
    for i = 1:j
        H(i,j) = w'*V(:,i);
        w = w - H(i,j)*V(:,i);
    end
    H(j+1,j) = norm(w);
    
    if abs(H(j+1,j)) < 1e-10
        sprintf('done without residual')        
        break;
    else
        V(:,j+1) = w./H(j+1,j);
    end
end
%接下来用qr分解求解最小二乘问题得到y：min||beta*e1-H*y||
e1=zeros(j+1,1);
e1(1)=beta;
[Q,R]=qr(H(1:j+1,1:j)); %R为m+1*m
R=R(1:j,1:j);
bb=Q'*e1;
y=backward(R,bb(1:j));
x=x0+V(:,1:j)*y;

r=norm(b-A*x);
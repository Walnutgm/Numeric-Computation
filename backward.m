function b=backward(A,b)
%回代法求解上三角形方程组
%输入A:一个上三角形方阵  b:右端向量
%输出：b为解
[n,~]=size(A);
for j=n:-1:2
    b(j)=b(j)/A(j,j);
    b(1:j-1)=b(1:j-1)-b(j)*A(1:j-1,j);
end
b(1)=b(1)/A(1,1);
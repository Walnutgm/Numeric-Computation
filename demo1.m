%使用共轭梯度法（CG）、预处理共轭梯度法（PCG）和广义最小残差法（GMRES）来求解A1x=b1
load('A1.mat')
load('b1.mat')
A=A1'*A1;
b=A1'*b1;
%由于矩阵A1不是对称矩阵，故在使用cg法和pcg法时采用求解与原方程组等价的方程组A'Ax=A'b
[x_cg,i_cg,r_cg]=CG(A,b,1e-4);

[x_pcg,i_pcg,r_pcg]=PCG(A,b,1e-4);

[x_gmres,r_gmres]=myGMRES(A1,b1); %运行时间较长

%使用共轭梯度法（CG）、预处理共轭梯度法（PCG）和广义最小残差法（GMRES）来求解A2x=b2
load('A2.mat')
load('b2.mat')
AA=A2'*A2;
bb=A2'*b2;
%由于矩阵A2不是对称矩阵，故在使用cg法和pcg法时采用求解与原方程组等价的方程组A'Ax=A'b
[x2_cg,i2_cg,r2_cg]=CG(AA,bb,1e-4);

[x2_pcg,i2_pcg,r2_pcg]=PCG(AA,bb,1e-4);

[x2_gmres,r2_gmres]=myGMRES(A2,b2);  %运行时间较长
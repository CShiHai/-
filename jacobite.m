function [k,Bk,V,D,Wc]=jacobite(A,jd,max1)
[n,n]=size(A);Vk=eye(n);Bk=A;state=1;k=0;P0=eye(n);
Aij=abs(Bk-diag(diag(Bk)));[m1 i]=max(Aij);
[m2 j]=max(m1);i=i(j);
while ((k<=max1)&(state==1))
    k=k+1,aij=abs(Bk-diag(diag(Bk)));[m1 i]=max(abs(aij));
    [m2 j]=max(m1);i=i(j),j,Aij=(Bk-diag(diag(Bk)));
    mk=m2*sign(Aij(i,j)),
    Wc=m2,Dk=diag(diag(Bk));Pk=P0;
    c=(Bk(j,j)-Bk(i,i))/(2*Bk(i,j)),
    t=sign(c)/(abs(c)+sqrt(1+c^2)),
    pii=1/( sqrt(1+t^2)), pij=t/( sqrt(1+t^2)),
    Pk(i,i)=pii;Pk(i,j)=pij;
    Pk(j,j)=pii; Pk(j,i)=-pij;
    Pk,B1=Pk'*Bk;B2=B1*Pk; Vk=Vk*Pk,Bk=B2,
if(Wc>jd)
    state=1;
else
       return
end
    Pk;Vk;Bk=B2;Wc;
end
if(k>max1)
    disp('请注意迭代次数k已经达到最大迭代次数max1,迭代次数k,对称矩阵Bk,以特征向量为列向量的矩阵V,特征值为对角元的对角矩阵D如下：')
else
    disp('请注意迭代次数k,对称矩阵Bk,以特征向量为列向量的矩阵V,特征值为对角元的对角矩阵D如下：')  
end
Wc;k=k; V=Vk;Bk=B2;D=diag(diag(Bk));
[V1,D1] =eig(A,'nobalance');

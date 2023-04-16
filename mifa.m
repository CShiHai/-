function [k,lambda,Vk,Wc]=mifa(A,V0,jd,max1)
lambda=0;k=1;Wc =1; jd=jd*0.1;state=1; V=V0;
while((k<=max1)&(state==1))
    Vk=A*V; [m,j]=max(abs(Vk)); mk=m; 
    tzw=abs(lambda-mk); Vk=(1/mk)*Vk;
    Txw=norm(V-Vk); Wc=max(Txw,tzw); V=Vk;lambda=mk;state=0;
    if(Wc>jd)
        state=1;
    end
    k=k+1;Wc=Wc;
end
if(Wc<=jd)
    disp('请注意：迭代次数k,主特征值的近似值lambda,主特征向量的近似向量Vk,相邻两次迭代的误差Wc如下：')  
else
    disp('请注意：迭代次数k已经达到最大迭代次数max1,主特征值的迭代值lambda,主特征向量的迭代向量Vk,相邻两次迭代的误差Wc如下：')     
end
 Vk=V;k=k-1;Wc;

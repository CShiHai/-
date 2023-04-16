function [x,h,k] = Jacobiiter(A,b,x0,N,p)  %N代表最大迭代次数
n=length(A);
n1=length(b);
D=zeros(n);L=zeros(n);U=zeros(n);
x=zeros(n,1);
x0=transpose(x0);
b=transpose(b);
if n~=n1
    disp('维数不一致')
    return
end
D(n,n)=A(n,n);
for i=1:n-1
    D(i,i)=A(i,i);
    for j=1:i-1
        L(i,j)=-A(i,j);
    end
    for r=i+1:n-1
        U(i,r)=-A(i,r);
    end
end
k=0;
while k<=N
    x=(eye(n)-inv(D)*A)*x0+inv(D)*b;
    k=k+1;
    h=norm((x-x0),inf);
    if h<=p
        break
    end
    
    x0=x;
    
end
if k>N
    disp(['迭代次数= ',num2str(k),',算法超过最大迭代次数。']);
    
else
    disp(['迭代次数= ',num2str(k)]);
    
end
x=x0
end
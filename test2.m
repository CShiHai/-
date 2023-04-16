clc
clear
A=[1 2 3;2 1 2;3 2 1];
b=[1;1;1];
n=size(A,1);%A的阶数
o=0.1;e=0.1;p=0.5;
[V,D] = eig(diag(b));
x=zeros(n,1);
y=zeros(n,1);
for i=1:n
    if D(i,i)~=0
        y(i,1)=1/D(i,i)/n;
    end
end
x=V*y;
for j=1:100
    w=zeros(n,1);
    s=0;
    I=[];
    a=A*x;
    if norm(x)>e
        for i=1:n
            w(i,1)=(p/2)*norm(x)^(p-2);
        end
    else
        for i=1:n
            w(i,1)=(p/2)*e^(p-2);
        end
    end
    Imin=min(o*w./b);
    umin=0-Imin;
    for i=1:n
        if o*w(i)/b(i)==Imin
            I=[I i];s=s+b(i)*a(i)/(umin*b(i)+o*w(i)).^2;
        end
    end
    for i=1:n
        if ismember(i,I)
            syms u f
            g1 = symsum(b(i)*a(i).^2,f,1,n);
            g = symsum(b(f)*a(f).^2/(u*b(f)+o*w(f)).^2,f,1,n)==1;
            answ = solve(g, u);
            x(i)=a(i)/(u*b(i)+o*w(i));
        else
            x(i)=a(i)/(umin*b(i)+o*w(i));
        end
    end
end


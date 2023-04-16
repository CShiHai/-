clc
clear
syms u
A1=[1 2 3;2 1 2;3 2 1];
b=[1;1;1];
% load model7;
% A1 = X;
% b = ones(200,1);
n=size(A1,1);%A的阶数
o=0.1;e=0.1;p=0.5;
B = diag(b);
[V,D] = eig(B);
[U,J] = eig(A1);
J = real(J);
U = real(U);
max1=J(1,1);s=1;
m = size(J,1);
optTol = 1e-9;
for i=2:m
    if J(i,i)>max1
        max1=J(i,i);
        s=i;
    end
end
s1 = U(:,s);
x=zeros(n,1);
y=zeros(n,1);
A = A1;
for i=1:n
    if D(i,i)~=0
        y(i,1)=1/D(i,i)/n;
    end
end
x=V*y;
for j=1:100
    xold = x;
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
            answ = solve(sum(b.*(a.*a)./((u.*b+o.*w).*(u.*b+o.*w)))==1, u);
            if size(answ,1)>1
                answmax=max(answ);
            else
                answmax=answ;
            end
            x(i)=a(i)/(answmax*b(i)+o*w(i));
        else
            x(i)=a(i)/(umin*b(i)+o*w(i));
        end
    end
    if norm(x-xold)<optTol
        CV= corr(x,s1);
        disp(x);
        disp(CV);
        break;
    end
end




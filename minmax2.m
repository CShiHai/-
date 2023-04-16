clc;
clear;

% A = [1 2 3;2 1 2;3 2 1];
load SpareMatrix1000; 
% A = X;


[n,m] = size(A);
e=0.1;p=0.1;
x=zeros(n,1);
for i = 1:n
    x(i)=1/(sqrt(n));
end
pe=p/(log(1+1/e));
for i = 1:10000
    w=1./(abs(x)+e);
    if pe < 2*max(abs(A*x)./w)
        B=abs(A*x)-pe/2*w;
        for j=1:n
            if B(j)<0
                B(j)=0;
            end
        end
        C=(abs(A*x)-pe/2*w).^2;
        for j=1:n
            if C(j)<0
                C(j)=0;
            end
        end
        x=(B.*sign(A*x))/sqrt(sum(C));
    else
        x=zeros(n,1);
    end
end

[U,J] = eig(A);
J = real(J);
U = real(U);
max1=J(1,1);s=1;
for i=2:m
    if J(i,i)>max1
        max1=J(i,i);
        s=i;
    end
end
s1 = U(:,s);


CV= corr(x,s1);
disp(abs(CV));
clc
clear
load SpareMatrix100;
A1=A;
% A1=[1 2 3;2 1 2;3 2 1];
n=size(A1,1);%A1的阶数
B=eye(n);
o=0.01;e=0.01;p=0.05;
[V,D] = eig(B);
[U,J] = eig(A1);
J = real(J);
U = real(U);
A = A1*J*A1'+B;
x=zeros(n,1);
y=zeros(n,1);
for i=1:n
    if D(i,i)~=0
        y(i,1)=1/sqrt(D(i,i))/n;
    end
end
x=V*y;
k=0;
w=zeros(n,1);


% options for convergence

optTol = 1e-9; % for convergence

for j=1:100
    if norm(x)>e
        for i=1:n
            w(i,1)=(p/2)*norm(x)^(p-2);
        end
    else
        for i=1:n
            w(i,1)=(p/2)*e^(p-2);
        end
    end
    [Q,W]=eig(A-o*diag(w));
    max=W(1,1);a=1;
    m = size(W,1);
    for i=2:m
        if W(i,i)>max
            max=W(i,i);a=i;
        end
    end
    max1=J(1,1);s=1;
    for i=2:m
        if J(i,i)>max1
            max1=J(i,i);
            s=i;
        end
    end
    s1 = U(:,s);
    xold = x;
    x=Q(:,a);
    if norm(x-xold)<optTol
        disp(j);
        CV= corr(x,s1);
        %CV = (q'*a1)/(norm(q)*norm(a1));
        disp(abs(CV));
%         c = x*x'-A;
%         p= norm(c,2);
%         disp(p);
        break;
    end
end
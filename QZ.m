clc
clear
% A=[1 2 3;2 1 2;3 2 1];
% b=[1;1;1];
tic;
load matrix500;
A = X;
n = size(A);
B=eye(n);
[V,D]=eig(A,B,'qz');
toc;
% disp(V);
% disp(D);

[U,J,W] = eig(A);
J = real(J);
U = real(U);
max1=J(1,1);
s=1;
m = n;
for i=2:m
    if J(i,i)>max1
        max1=J(i,i);
        s=i;
    end
end
s1 = U(:,s);

d=1;
for i=2:n
    if D(i,i)>max1
        max1=D(i,i);
        d=i;
    end
end
x = V(:,d);


CV= corr(x,s1);
%CV = (q'*a1)/(norm(q)*norm(a1));
disp(abs(CV));



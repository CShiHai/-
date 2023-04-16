clc;
clear;
tic;
% A=[1 2 3;2 1 2;3 2 1];
load matrix10;
A = X;
n = size(A);
B=eye(n);
itermax = 1000;
[d,v,history,historyend,numrot]=jacobi(A,itermax);

a=1;
max1=d(1);
for i=2:n
    if d(i)>max1
        max1=d(i);
        a=i;
    end
end
x = v(:,a);
toc;

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


CV= corr(x,s1);
%CV = (q'*a1)/(norm(q)*norm(a1));
disp(abs(CV));
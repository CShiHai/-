% A = [-6.9,14,0;
%      -5,10.1,0;
%      -1,0,-0.1];
clc;
clear;
tic;
% A=[1 2 3;2 1 2;3 2 1];
% b=[1;1;1];
% B = diag(b);
load matrix600;
A = X;
n = size(A);
B=eye(n);
A = B^-1*A;
N=100;
ep=1e-4;
n=length(A);
y=ones(n,1);
k=0;
m1=0;
while k<=N
   x=A*y;
   m=max(abs(x));
   y=x/m;
   if abs(m-m1)<ep
        break;
   end          
   m1=m;
   k=k+1;
end

final_answer_u = m;
final_answer_x = x;
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

CV= corr(final_answer_x,s1);
%CV = (q'*a1)/(norm(q)*norm(a1));
disp(abs(CV));


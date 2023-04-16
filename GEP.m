clc
clear
A=[1 2 3;2 1 2;3 2 1];
b=[1;1;1];
B = diag(b);
[V,D,W]=eig(A,B);
disp(V);
disp(D)

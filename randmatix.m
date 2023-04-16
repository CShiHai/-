clc;
clear;
m = 500; % number of rows
n = 500; % number of columns
C = randn(m, m);
A = C'*C;
X =A; 
save('matrix600','X')

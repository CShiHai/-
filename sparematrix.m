clc;
clear;
m = 100; % number of rows in the first 10 elements
n = 500; % number of columns
A = zeros(n, n); % generate a matrix with m+10 rows

% set the first m rows of each column to follow normal distribution
A(1:m, :) = randn(m, n);

% set the remaining rows to 0
A(m+1:end, :) = 0;
save('SpareMatrix502','A')
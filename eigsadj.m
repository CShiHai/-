% adjusted generalized eigenvalue decomposition 

function [T1,D1] = eigsadj(A1,B1,k)

if nargin == 2
[T1, D1] = eigs(A1, B1);
else  
[T1, D1] = eigs(A1, B1,k);  
end

tmp = diag(T1'*B1*T1);
for kk = 1:k
    T1(:,kk) = T1(:,kk) ./ sqrt(tmp(kk));
end
end
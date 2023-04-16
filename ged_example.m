%% Lobular Freeze data analysis example (PCA)
clear;clc;
load SpareMatrix1000; 
% X = [1 2 3;2 1 2;3 2 1];
X = A;
[n,nvariables] = size(X); 
p = n; 
rng(1);
trainsize = floor( n / 3 );
ID = zeros(n,1);[~,rind]=sort(rand(n,1));
itrain = logical(ID); itrain(rind(1:trainsize)) = true;

X = X(:,1:p); 
V1 = X;
X = [X randn(n,p)*2]; 
[n,p]=size(X);
X = X - repmat(mean(X),size(X,1),1);
X1 = X(itrain,:);
An = cov(X1);

d = n; 

[U,J] = eig(An);max1=J(1,1);s=1;
m = size(An,1);
for i=2:m
    if J(i,i)>max1
       max1=J(i,i);
       s=i;
    end
end
 s1 = U(:,s);
 
[B,A,k,eigenvectors,eigenvalues,lambda] = seig(V1,An,[],d);
CV= corr(s1,eigenvectors(:,1));
% CV = (eigenvectors(:,1)'*s1)/(norm(eigenvectors(:,1))*norm(s1));
disp(abs(CV));
% [CV,lvec,Vcell, Vonly] = POIcv(B, A, An,B, k);
% disp(CV);
% disp(lvec);
% disp(Vcell);
% disp(Vonly);


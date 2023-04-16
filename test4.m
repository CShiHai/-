clc
clear
load model7;
A1=X;
n = size(A1,1);
x = zeros(n,1);
o=0.1;e=0.1;p=0.5;
pw = p/log(1+e.^-1);
for l=1:n
    for i=0:n-1
        w(i)=(x(i)+e).^-1;
    end
    a
    
end


clc
clear
A1=[1 2 3;2 1 2;3 2 1];
% load model7;
% A1 = X;
[U,J] = eig(A1);
J = real(J);
U = real(U);
n = size(A1,1);
x = zeros(n,1);
w = zeros(n,1);
o=0.1;e=0.1;p=0.5;
pw = p/log(1+e.^-1);
max1=J(1,1);s=1;
for i=2:n
    if J(i,i)>max1
        max1=J(i,i);
        s=i;
    end
end
s1 = U(:,s);
xold = x;
optTol = 1e-9;
for i=1:n
    w(i)=(abs(x(i))+e).^-1;
end
I=abs(A1*x).*(w).^(-1);
Imax=max(abs(A1*x).*(w).^(-1));
[row,col]=find (I==Imax);
rowmin=max(row);
if pw<2*rowmin
    B=A1*x;
    for i=1:n
        %             x(i)=((abs(B(i))-(pw/2).*w(i))+sign(B(i)))/sqrt(sum((abs(B)-(pw/2).*w(i))).^2);
        x(i)=((abs(B(i))-(pw/2).*w(i))+sign(B(i)))/sqrt(sum((abs(B)-(pw/2).*w(i)).^2));
    end
else
    x = zeros(n,1);
end
while norm(x-xold)>optTol
    for i=1:n
        w(i)=(abs(x(i))+e).^-1;
    end
    I=abs(A1*x).*(w).^(-1);
    Imax=max(abs(A1*x).*(w).^(-1));
    [row,col]=find (I==Imax);
    rowmin=max(row);
    if pw<2*rowmin
        B=A1*x;
        for i=1:n
            %             x(i)=((abs(B(i))-(pw/2).*w(i))+sign(B(i)))/sqrt(sum((abs(B)-(pw/2).*w(i))).^2);
            x(i)=((abs(B(i))-(pw/2).*w(i))+sign(B(i)))/sqrt(sum((abs(B)-(pw/2).*w(i)).^2));
        end
    else
        x = zeros(n,1);
    end
end
CV= corr(x,s1);
disp(CV);



%An highlighted block
function x=jacobi_fun(A,b,n,x0,tol,N)
x=zeros(n,1);      % 给x赋值
k=0;
while k<N
    for i=1:n
     x(i)=(b(i)-A(i,[1:i-1,i+1:n])*x0([1:i-1,i+1:n]))/A(i,i);
    end
    if norm(x-x0)<tol
        break;
    end
    x0=x;
    k=k+1;

    disp(['when k=',num2str(k)])
    disp('x=');
    disp(x);                       %输出计算的中间结果
end

if k==N
    disp('迭代次数已到达上限!');
end
disp(['迭代次数 k=',num2str(k)])

end

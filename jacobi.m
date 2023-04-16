function [d,v,history,historyend,numrot]=jacobi(a_in,itermax)
% [d,v,history,historyend,numrot]=jacobi(a_in,itermax)
% computes the eigenvalues d and
% eigenvectors v of the real symmetric matrix a_in,
% using rutishausers modfications of the classical
% jacobi rotation method with treshold pivoting.
% history(1:historyend) is a column vector of the length of
% total sweeps used containing the sum of squares of
% strict upper diagonal elements of a. a is not
% touched but copied locally
% the upper triangle is used only
% itermax is the maximum number of total sweeps allowed
% numrot is the number of rotations applied in total

% check arguments
siz=size(a_in);
if siz(1) ~= siz(2)
    error('jacobi : matrix must be square ' );
end
if norm(a_in-a_in',inf) ~= 0
    error('jacobi ; matrix must be symmetric ');
end
if ~isreal(a_in)
    error(' jacobi : valid for real matrices only');
end
if nargin==1
    itermax=20;
end
    
n=siz(1);
v=eye(n);
a=a_in;
history=zeros(itermax,1);
d=diag(a);
bw=d;
zw=zeros(n,1);
iter=0;
numrot=0;
while iter < itermax
iter=iter+1;
history(iter)=sqrt(sum(sum(triu(a,1).^2)));
historyend=iter;
tresh=history(iter)/(4*n);
if tresh ==0
return;
end
for p=1:n
for q=p+1:n
gapq=10*abs(a(p,q));
termp=gapq+abs(d(p));
termq=gapq+abs(d(q));
if iter>4 && termp==abs(d(p)) && termq==abs(d(q))
% annihilate tiny elements
a(p,q)=0;
else
if abs(a(p,q)) >= tresh
%apply rotation
h=d(q)-d(p);
term=abs(h)+gapq;
if term == abs(h)
t=a(p,q)/h;
else
theta=0.5*h/a(p,q);
t=1/(abs(theta)+sqrt(1+theta^2));
if theta < 0
t=-t;
end
end
c=1/sqrt(1+t^2);
s=t*c;
tau=s/(1+c);
h=t*a(p,q);
zw(p)=zw(p)-h; %accumulate corrections to diagonal elements
zw(q)=zw(q)+h;
d(p)=d(p)-h;
d(q)=d(q)+h;
a(p,q)=0;
%rotate, use information from the upper triangle of a only
%for a pipelined cpu it may be better to work
%on full rows and columns instead
for j=1:p-1
g=a(j,p);
h=a(j,q);
a(j,p)=g-s*(h+g*tau);
a(j,q)=h+s*(g-h*tau);
end
for j=p+1:q-1
g=a(p,j);
h=a(j,q);
a(p,j)=g-s*(h+g*tau);
a(j,q)=h+s*(g-h*tau);
end
for j=q+1:n
g=a(p,j);
h=a(q,j);
a(p,j)=g-s*(h+g*tau);
a(q,j)=h+s*(g-h*tau);
end
% accumulate information in eigenvectormatrix
for j=1:n
g=v(j,p);
h=v(j,q);
v(j,p)=g-s*(h+g*tau);
v(j,q)=h+s*(g-h*tau);
end
numrot=numrot+1;
end
end %if
end % for q
end % for p
bw=bw+zw;
d=bw;
zw=zeros(n,1);
end %while



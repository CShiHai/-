function [B,A,k,v,d,lambda,q] = seig(V1,A,B,k,sp,option,abssp)
% SEIG Sparse Generalized eigendecomposition via Penalized Orthogonal
% Iteration (POI)
%
% This is a wrapper for a simple use of POI algorithm to solve a GEP. 
%这是简单使用POI算法求解GEP的包装
% 
% [v,d,lambda,q] = seig(A,B) ;
% [v,d,lambda,q] = seig(A,[]);
% [v,d,lambda,q] = seig(A);   A and B are symmetric non-negative definite matrices,
%  A和B是对称的非负定矩阵,大小为p*p,v是第一个特征向量，d是对应的最大的特征值。
% lamba 是使用的调整参数
% of the size p x p, v is the first eigenvector, d is the corresponding,
% largest eigenvalue, lambda is the tuning parameter used.
% If B is the identity matrix, eye(p), or is empty,[], then this corresponds
% to the penalized solution to the eigenvalue decomposition. Here, q is the 
% orthogonal basis that spans columns of v (Note that v is only B-orthogonal).
%
%如果B是单位矩阵eye（p）或为空[]，则这对应于特征值分解的惩罚解。
%这里，q是跨越v列的正交基（注意v只是B正交)
%
% [v,d,lambda,q] = seig(A,B,k); k is the dimension of the subspace, or the number of
% k是子空间的维数或特征向量的数量。默认k=1，对于k>1,
% v是大小p x k的特征向量矩阵估计，d是特征值矩阵，k x k对角矩阵由递减特征值估计组成。
% eigenvectors. By default, k = 1. For k > 1, v is the eigenvector matrix
% estimate of size p x k, and d is the eigenvalue matrix, the k x k diagonal
% matrix consisting of decreasing eigenvalue estimates. Here, q is the 
% orthogonal basis that spans columns of v (Note that v is only B-orthogonal).
% 
% [v,d] = seig(A,B,k,sp); sp is a "scaled" tuning parameter used in the
% penalization of POI. (sp stands for a sparsity parameter.) 
%sp是用于惩罚POI的“缩放”调整参数。（sp代表稀疏度参数。）
% For sp = 0, the tuning parameter lambda = 0. For sp in (0,1], lambda = sp
% * lambda_max, where lambda_max is the maximum nontrivial value of lambda, 
% obtained by POIlim.m. By default, sp = 1 / 2. 
%对于sp=0，调谐参数lambda=0。对于sp in（0,1]），lambda=sp*lambda_max，
% 其中lambda_max是由POIlim.m获得的lambda的最大非平凡值。默认情况下，sp=1/2。

% If sp = [], then the default value is used. We strongly 
% recommend to use the cross-validation in POIcv.m to tune lambda. 
%如果sp=[]，则使用默认值。我们强烈建议使用POIcv.m中的交叉验证来调整lambda
% [v,d] = seig(A,B,k,sp,option), where option is a string designating
% the type of POI algorithm. By default, option = 'POI-C'. 
%其中option是指定POI算法类型的字符串。默认情况下，选项=“POI-C”。
%  option =  'POI-L', 'POI-C', 'FastPOI-L', or 'FastPOI-C';
%   A = Penalized Orthogonal Iteration with coefficient-wise L1 penalty (lasso)
%   C = Penalized Orthogonal Iteration with coordinate-wise L1 penalty (group lasso)
%   Da= Fast POI with coefficient-wise L1 penalty (lasso)
%   D =  Fast POI with coordinate-wise L1 penalty (group lasso)
%  
% See also POI, POIlim, POIv, POIcv. 
%
% Last updated May 2018
% Sungkyu Jung

default.option = 'POI-C';
default.k = 1; 
% if nargin < 2 
%     B = [];
% end
% if nargin < 3 
%     k = default.k; 
%     option = default.option;
%     sp = 1/2 ; 
% elseif nargin < 4
%     option = default.option;
%     sp = 1/2 ; 
% elseif nargin < 5
%     option = default.option;
% end

if nargin < 3 
    B = [];
end
if nargin < 4 
    k = default.k; 
    option = default.option;
    sp = 1/2 ; 
elseif nargin < 5
    option = default.option;
    sp = 1/2 ; 
elseif nargin < 6
    option = default.option;
end

if isempty(sp)
    sp = 1/2 ; 
end
lambda = POIlim(A,option,k)*sp; 


if nargin == 6  % debug purpose only 
    if abssp; lambda = sp; end
end


if isempty(B)
    B = eye(size(A,1));
else
    rB = rank(B);
    if rB < size(B,1)
        p = size(B,1); 
        eps = log(p) / rB ; 
        eps2 = min(diag(B)) ;
        if  eps2 > 0 
            eps = min(eps,eps2);
        end
        B = B + eps * eye(p);
    end
end




maxIterInner =500;
Q = POI(B, A, lambda, k, option, maxIterInner);
%Q = POI(B, A, lambda, k);
gepsolutions=POIv(B,A,Q);
%gepsolutions=POI(B,A,100,100,POI-L,A,100,[]);
v = gepsolutions.U;
% disp(v);
% a = Q*Q'-A;
% p= norm(a,2);
% disp(p);
d = gepsolutions.Lambda;
q = gepsolutions.Q;
% inputstruct.nGrid=11;
% inputstruct.maxIterInner=100;
% inputstruct.lvec=lambda;


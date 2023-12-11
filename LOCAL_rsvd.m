function [U,D,V] = LOCAL_rsvd(A,k,p)

n         = size(A,2);
ell       = k + p;
Omega     = randn(n,ell);
Y         = A*Omega;
[Q,~,~]   = qr(Y,0);
B         = Q'*A;
[UU,D,V]  = svd(B,'econ');
U         = Q*UU(:,1:k);
D         = D(1:k,1:k);
V         = V(:,1:k);

return
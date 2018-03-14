function [ S ] = fcls( A, X )

if (~ismatrix(A)), error('A must be a L by M matrix.'); end
if (~ismatrix(X)), error('X must be a L by N matrix.'); end

[L1,M] = size(A);
[L2,N] = size(X);

if (L1 ~= L2), error('The number of bands is not the same.'); end

S = zeros(M, N);
A_ = [1e-5*A;ones(1,M)];
X_ = [1e-5*X;ones(1,N)];
for i=1:N, S(:,i) = lsqnonneg(A_,X_(:,i)); end
S = min(max(S,eps),1.0000000);

end


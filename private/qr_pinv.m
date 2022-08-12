function Xpinv = qr_pinv(X, eps)
% Efficiently compute the Moore-Penrose pseudoinverse of X using QR
% decomposition. Faster than Matlab's built-in pinv(), which uses
% singular value decomposition. However, qr_pinv() requires that X
% have full column rank.
%
% X must be full column rank to precision eps (default eps = 0.0001).

% Compute the QR decomposition of X.
[Q, R] = qr(X, 0);

% Check that X is full rank.
if nargin < 2
    eps = 0.0001; % default precision for determining rank
end
if size(X,2) > size(X,1) || any(abs(diag(R)) < eps)
    error('X is not full rank.');
end

% Obtain the pseudoinverse of X by solving the following system using
% back-substitution of the upper triangular matrix R:
% R * Xpinv = Q'
% Note that we only need the non-singular parts of Q and R.
Xpinv = linsolve(R(1:size(X,2),1:size(X,2)), Q(:,1:size(X,2))', struct('UT', true));

end

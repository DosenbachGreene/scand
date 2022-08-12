function SE = covB_to_SE(covB)
% Compute standard error from variance-covariance matrix.
%
% Helper function to compute the standard errors of fixed effect estimates
% (Betas, B) from the variance-covariance matrix of B, covB.  The full covB
% is required when performing a Wald test. The standard errors are
% sufficient for performing a t-test. The standard errors are computed as:
%
% SE = sqrt(diag(covB))
%
% SE has one row for each row in B or column in the design matrix X. When
% covB is a 3-dimensional matrix such that the third dimension corresponds
% to independent columns in the matrix of observations Y then SE will have
% one column for each column in Y.
%
% See also: SCAND

% Make sure argument is a square matrix.
assert(size(covB,1) == size(covB,2));

if ndims(covB) < 3 || size(covB,3) == 1
    SE = sqrt(diag(covB));
else
    SE = zeros(size(covB,1), size(covB,3));
    parfor i=1:size(covB,3)
        SE(:,i) = sqrt(diag(covB(:,:,i)));
    end
end

end
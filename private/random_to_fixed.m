function Xaug = random_to_fixed(Z, G)
% Generate fixed effects covariates to model random effects.
%
% Meant to be used as part of SCAND.  Takes the random effects design
% matrix Z and grouping G passed to SCAND.  Generates a set of fixed
% effect covariates to model the random effects.  The output, Xaug, is
% meant to augment the fixed effects design matrix X.  The columns of Xaug
% are centered such that the fixed effects in X model the "overall"
% (marginal) effects after accounting for group/cluster membership.
%
% Note: Does not guarantee degrees of freedom.  The caller should verify:
% assert(size(X,2) + size(Xaug,2) < size(X,1))

% Sanity check arguments.
if ~iscell(Z) && ~iscell(G)
    Z = {Z};
    G = {G};
end
if iscell(Z) && iscell(G)
    assert(length(Z) == length(G))
    for i=1:length(Z)
        assert(size(Z{i},2) < size(Z{i},1));
        assert(size(G{i},1) == size(Z{i},1));
        assert(iscategorical(G{i}));
        assert(size(G{i},2) == 1);
    end
else
    error('If one of Z or G is a cell array then both must be cell arrays.');
end

% Allocate a matrix for augmenting X with columns to control for the
% random effects.
Xaug_ncols = 0; % number of columns to augment with
for i_GZ=1:length(G)
   Xaug_ncols = Xaug_ncols + size(Z{i_GZ},2) * (length(categories(G{i_GZ})) - 1);
end
Xaug = zeros(size(Z{1},1), Xaug_ncols);

% Generate fixed effects predictors for the "random" effects.
% Each column in Z needs length(cats)-1 predictors where cats is the
% number of categories in the corresponding group G. The predictors are
% centered so that they will not affect the estimates of the fixed
% effects except in the (desired) case where the fixed and random
% effects are correlated.
%
% For example, suppose there is a random effect with two categorical
% levels.  We would augment X with 2-1 = 1 "random" predcitor. The
% predictor would look like [1; 1; 1; ... -1; -1; -1] where 1 encodes
% the first category and -1 encodes the second category.
%
% For a random effect with three categorical levels we would augment X
% with 3-1 = 2 "random" predictors.  The predictors would look like:
%    1    -0.5
%    1    -0.5
%   -0.5   1
%   -0.5   1
%   -0.5  -0.5
%   -0.5  -0.5
% Where [1, -0.5] encodes the first category, [-0.5, 1] encodes the
% second category, and [-0.5, -0.5] encodes the third category. In
% general, if there are n categories then there are n-1 predictors
% that take values of 1 or -1/(n-1).
Xaug_col = 1;
for i_GZ=1:length(G) % pairs of groups and random effect matrices
    cats = categories(G{i_GZ});
    c = 1/(length(cats)-1); % value for encoding categories
    for j_cat=1:length(cats) % categories within groups
        for k_Z=1:size(Z{i_GZ},2) % columns within random effect matrices
            % Need length(cats)-1 predictors, so skip last category.
            if j_cat < length(cats)
                % Generate a predictor column.
                Xaug(:,Xaug_col) = Z{i_GZ}(:,k_Z) .* ((G{i_GZ} == cats(j_cat)) .* (1+c) - c);
                Xaug_col = Xaug_col + 1;
            end
        end
    end
end

end
function covB = swe_block(Xpinv, resid, G, SE_only)
% Covariance matrix of the fixed effects using the sandwich estimator.
%
% Compute the variance-covariance matrix of some marginal, fixed effects B
% using the Huber-White Sandwich Estimator under the assumption that the
% variance-covariance matrix V of the residuals has a block-diagonal
% structure.
%
% Briefly, the sanwich estimator is given by:
% B = inv(X'*W*X)*X'*W*Y
% covB = inv(X'*W*X) * X'*W*V*W*X * inv(X'*W*X)
% Where V is the variance-covariance matrix of the residuals and and W is a
% weighting matrix.  In a mixed effects model the weighting matrix
% W = inv(V). However, computing V and its inverse can be computationally
% expensive. It turns out that the sandwich estimator produces
% asymptotically-correct estimates of covB for empirical estimates of V
% even when W is misspecified. Therefore we take:
% W = I
% V = resid * resid'
% Which leads to:
% covB = Xpinv*V*Xpinv', Xpinv = pinv(X)
%
% However, the above formulation leads to covB = 0 unless the structure of
% V is constrained. In this block version of the SwE we constrain V to be
% zero between any two residuals that do not share a group in G.
%
% For further details, refer to the software package:
% https://www.nisox.org/Software/SwE/
% And to Bryan Guillaume thesis work:
% Guillaume B, Hua X, Thompson PM, Waldorp L, Nichols TE; Alzheimer's Disease Neuroimaging Initiative. Fast and accurate modelling of longitudinal and repeated measures neuroimaging data. Neuroimage. 2014 Jul 1;94:287-302. doi: 10.1016/j.neuroimage.2014.03.029
%
% Also see commentary by David Freedman for caveats:
% Freedman DA. On The So-Called "Huber Sandwich Estimator" and "Robust Standard Errors." The American Statistician. 2006;60(4):299-302. doi:10.1198/000313006x152207

% Perform sanity checks.
assert(size(Xpinv,2) == size(resid,1));
assert(size(Xpinv,2) > size(Xpinv,1));
if ~iscell(G)
    G = {G};
end
for i=1:length(G)
    assert(size(G{i},1) == size(Xpinv,2));
    assert(iscategorical(G{i}));
    assert(size(G{i},2) == 1);
end
if nargin == 3
    SE_only = false;
end
assert(islogical(SE_only));

% Build up the structure of V.
% V = false(size(Xpinv,2));
% for i_G=1:length(G)
%     cats = categories(G{i_G});
%     Z = false(size(Xpinv,2), length(cats));
%     for j_cat=1:length(cats)
%         Z(G{i_G} == cats(j_cat), j_cat) = true;
%     end
%     V = V | (Z*Z');
% end

% Find combinations of categories.

% GG = G{1};
% for i_G=2:length(G)
%     GG = GG .* G{i_G};
% end
% cats = categories(GG);

%% Find independent blocks.

% Allocate a vector of integers, one for each row in G.
% A value of zero indicates the row is not yet assigned to a block.
blocks = zeros(length(G{1}),1, 'int32');

% The number of independent blocks.
nblocks = int32(0);

% Iterate over the rows in G and assign each row to a block.
for i=1:length(G{1})
    % Skip rows that have already been assigned.
    if blocks(i) == 0
        % Increment the number of blocks.
        nblocks = nblocks + 1;
        
        % Assign this row to this block.
        blocks(i) = nblocks;
        
        % Assign all matching rows in G to this block.
        % Repeat until there are no more new rows to assign to this block.
        try_assign_more = true;
        while try_assign_more
            % Default to not assigning more rows to the block on the next
            % iteration.
            try_assign_more = false;
            
            % Iterate over columns in G.
            for j=1:length(G)
                % Find rows in this column of G with the same value as any
                % of those rows in this column of G which are already in
                % the block.
                old_block = blocks == nblocks;
                new_block = ismember(G{j}, G{j}(blocks == nblocks));
                
                % Have we found new rows to add to the block?
                if sum(new_block) > sum(old_block)
                    % Add the new rows to the block.
                    blocks(new_block) = nblocks;
                    % Try to find more rows to add to the block until there
                    % are no new rows added.
                    try_assign_more = true;
                end
            end
        end
    end
end

if ndims(resid) == 1 || size(resid,2) == 1
    % Allocate memory.
    if SE_only
        covB = zeros(size(Xpinv,1), 1);
    else
        covB = zeros(size(Xpinv,1), size(Xpinv,1));
    end
    
    % We could just take covB = Xpinv * (resid * resid') * Xpinv'
    % For performance, we will instead compute the covB for each
    % independent block and sum them together.
    for i=1:nblocks
        this_block = blocks == i;
        half_sandwich = Xpinv(:,this_block)*resid(this_block);
        if SE_only
            % We only need the diagonal of covB, which we can compute more
            % efficiently as the column sum of element-wise multiplication.
            % Equivalent to:
            % covB = covB + diag(half_sandwich * half_sandwich'
            covB = covB + sum(half_sandwich .* half_sandwich, 2);
        else
            % Compute the sandwich estimator for this block.
            covB = covB + half_sandwich * half_sandwich';
        end
    end
else
    % Compute a covariance matrix for each column of residuals.
    if SE_only
        covB = zeros(size(Xpinv,1), size(resid,2));
    else
        covB = zeros(size(Xpinv,1), size(Xpinv,1), size(resid,2));
    end
    if SE_only
        parfor i_resid=1:size(resid,2)
            for j_block=1:nblocks
                this_block = blocks == j_block;
                half_sandwich = Xpinv(:,this_block)*resid(this_block,i_resid);
                covB(:,i_resid) = covB(:,i_resid) + sum(half_sandwich .* half_sandwich, 2);
            end
            if mod(i_resid, 1000) == 0
                fprintf(1, '%0.2f%%\n', i_resid/size(resid,2)*100);
            end
        end
    else
        parfor i_resid=1:size(resid,2)
            for j_block=1:nblocks
                this_block = blocks == j_block;
                half_sandwich = Xpinv(:,this_block)*resid(this_block,i_resid);
                covB(:,:,i_resid) = covB(:,:,i_resid) + half_sandwich * half_sandwich';
            end
            if mod(i_resid, 10000) == 0
                fprintf(1, '%0.2f%%\n', i_resid/size(resid,2)*100);
            end
        end
    end
end

% If we are returning standard errors instead of covariance then we need to
% take the square root.
if SE_only
    covB = sqrt(covB);
end

end
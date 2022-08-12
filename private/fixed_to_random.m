function u = fixed_to_random(Z, G, Baug)
% Generate random effects from fixed estimates of random covariates.
%
% Meant to be used as part of SCAND. Takes the random effects design matrix
% Z and grouping G passed to SCAND and the estimates (columns in B)
% corresponding to Xaug returned by random_to_fixed(). The output, u, are
% approximations of the random effects such that the cluster-dependent
% estimates are B + u. The u are returned in the same order as
% FITLMEMATRIX.  The random effects are ordered first by G, then by
% categories within G, and lastly by columns in Z.
%
% If Baug has multiple columns for independent columns in Y then u will
% have the same number of columns.

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

% Allocate memory for random effects.
u_nrows = 0;
for i_GZ=1:length(G)
   u_nrows = u_nrows + size(Z{i_GZ},2)*length(categories(G{i_GZ})); 
end
u = zeros(u_nrows, size(Baug, 2));

Bcol = 1; % start at first column in Baug
urow = 1;
for i_GZ=1:length(G) % pairs of groups and random effect matrices
    cats = categories(G{i_GZ});
    c = 1/(length(cats)-1); % same c used when encoding categories above
    % Number of random columns in B corresponding to this grouping.
    nBcols = size(Z{i_GZ},2) * (length(cats)-1);
    for j_cat=1:length(cats) % categories within groups
        for k_Z=1:size(Z{i_GZ},2) % columns within random effect
            % Columns in B corresponding to this column in Z.
            Bcols = (Bcol + k_Z - 1) : size(Z{i_GZ},2) : (Bcol + nBcols - 1);
            if j_cat < length(cats)
                % The random effect is:
                % u = -c*B1 + -c*B2 + ... + Bz + -c*B3  + -c*B4 + ...
                % Where Bz is for *this* column in Z and B1... are
                % for the other columns in Z for this group (i.e.
                % the columns for the other categories).
                u(urow, :) = sum(Baug(Bcols([1:(j_cat-1), (j_cat+1):end]),:)) * -c + Baug(Bcols(j_cat),:);
            else
                % Treat the last category differently.
                % The random effect is just the sum of all the
                % corresponding columns in Z times the
                % rescaling/indexing factor:
                % u = -c*B1 + -c*B2 + ...
                u(urow, :) = sum(Baug(Bcols,:)) * -c;
            end
            urow = urow + 1;
        end
    end
    Bcol = Bcol + nBcols;
end

end
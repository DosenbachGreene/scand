function model = scand(X, Y, Z, G, varargin)
% SCAND Approximately fit a linear mixed effects model to the data.
%
% SCAND (Speedy Cluster Analysis for Neuroimaging Data) is a high-level
% function for approximating a linear mixed effects model to large data
% for which optimal fitting with REML (restricted maximum likelihood)
% algorithms would be prohibitively slow. The function syntax approximates
% that of Matlab's FITLMEMATRIX:
%
%   model = scand(X, Y, Z, G, 'method', method_name ...)
%
% And approximately fits a subset of the random effects models supported by
% fitlmematrix, specifically:
%
%   Y = X*B + sum(sigma2_i * Z_i * Z_i') + e
%
% The required parameters are:
%
%   X: The n observations x p predictors fixed-effects design matrix.
%   Y: The n x m observed data, with m independent columns.
%   Z: The n x q predictors random-effects design matrix, see notes below.
%   G: An n x 1 categorical grouping variable.
%
% Arguments X and Y have the same meaning as in ordinary least squares
% regression. The random-effects design matrix Z follows a Matlab
% convention and is different from the "raw" fixed effects design matrix.
% Consider the simple case of modeling a random interecept for three
% clusters (groups). A typical random-effects design matrix strcmpwould have
% three columns encoding the three groups with ones and zeros:
%
%   Z =
%        1     0     0
%        1     0     0
%        0     1     0
%        0     1     0
%        0     0     1
%        0     0     1
%
% Matlab (and scand) use an alternative convention in which Z is a single
% column and groups are encoded in a separate, categorical variable G:
%
%   Z =
%        1
%        1
%        1
%        1
%        1
%        1
%   G = 
%        1 
%        1 
%        2 
%        2 
%        3 
%        3 
%
% If there is more than one grouping variable, then Z and G can each be a
% cell array of the same length, where each cell contains a matrix Z and
% vector G.
% 
% SCAND accepts an optional argument, method, which can be one of:
%
%   swe
%       Fits random effects as fixed effects and adjusts the standard error
%       of the marginal model using a sandwich estimator. (default)
%   
%   mom
%       Method of moments, slower than swe but also generates estimates of
%       the random effect variances. May also generate more efficient
%       estimates of the marginal standard errors under some circumstances.
%
%   reml
%       Calls Matlab's FITLMEMATRIX. Slowest but most accurate. Useful for
%       comparing models fit with one of the approximate methods above
%       against the "gold standard" results. Any additional arguments
%       passed to SCAND will be forwarded to FITLMEMATRIX.
%
% The output of SCAND is a structure with the following elements:
%
%   B
%       The (marginal) fixed effects estimates,
%       in the same order as the columns in X.
%
%   covB
%       The variance-covariance matrix of the fixed effects, B.
%       Note: to compute t-values, take B ./ sqrt(diag(covB))
%
%   SE
%       If SE_only is true (see below), returns only the standard error
%       instead of covB.  SE = sqrt(diag(covB))
%
%   u
%       The random effects estimates in the same order returned by Matlab's
%       randomEffects() function.  The random effects are ordered first by
%       G, then by categories within G, and lastly by columns in Z.
%
% For the method swe, SCAND accepts another optional argument, SE_only
% (default = false).  When SE_only is true then a more performant algorithm
% is used to compute only the standard errors of covB: SE =
% sqrt(diag(covB)).  This is useful for performing t-tests vs Wald tests.
%
% See also: covB_to_SE()

% Parse arguments.
p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'X');
addRequired(p, 'Y');
addRequired(p, 'Z');
addRequired(p, 'G');
addParameter(p, 'method', 'swe');
addParameter(p, 'SE_only', false);
parse(p, X, Y, Z, G, varargin{:});
X = p.Results.X;
assert(size(X,2) < size(X,1));
Y = p.Results.Y;
assert(size(Y,1) == size(X,1));
Z = p.Results.Z;
G = p.Results.G;
if ~iscell(Z) && ~iscell(G)
    Z = {Z};
    G = {G};
end
if iscell(Z) && iscell(G)
    assert(length(Z) == length(G))
    for i=1:length(Z)
        assert(size(Z{i},1) == size(X,1));
        assert(size(Z{i},2) < size(Z{i},1));
        assert(size(G{i},1) == size(X,1));
        assert(iscategorical(G{i}));
        assert(size(G{i},2) == 1);
    end
else
    error('If one of Z or G is a cell array then both must be cell arrays.');
end
method = p.Results.method;
assert( strcmp(method, 'swe') | strcmp(method, 'mom') | strcmp(method, 'reml') );
reml_args = p.Unmatched;
if ~strcmp(method, 'reml')
    assert(isempty(fieldnames(reml_args)));
end
SE_only = p.Results.SE_only;
assert(islogical(SE_only));
if SE_only
    assert(strcmp(method, 'swe'));
end
reml_args = struct2opt(reml_args);
clear p

% Invoke the selected method.
if strcmp(method, 'reml')
    % Easy case, just delegate to fitlmematrix.
    if ndims(Y) == 1 || size(Y,2) == 1
        % Just one REML model to fit.
        model = fitlmematrix(X, Y, Z, G, reml_args{:});
        % Fixed effects.
        B = fixedEffects(model);
        % Covariance of fixed effects.
        covB = model.CoefficientCovariance;
        % Random effects.
        u = randomEffects(model);
        % Variance of each random effect in Z.
        % The mse is the remaining homoskedastic variance.
        [psi, mse] = covarianceParameters(model);
        sigmas = cell2mat(cellfun(@(x)diag(x), psi, 'UniformOutput', false));
    else
        % Preallocate memory for results.
        B = zeros(size(X,2), size(Y,2));
        covB = zeros(size(X,2), size(X,2), size(Y,2));
        u_len = 0; % number of random effects to allocate memory for
        z_len = 0; % number of random effect variances (sigmas)
        for i=1:length(G)
           u_len = u_len + size(Z{i},2)*length(categories(G{i})); 
           z_len = z_len + size(Z{i},2);
        end
        u = zeros(u_len, size(Y,2));
        sigmas = zeros(z_len, size(Y,2));
        mse = zeros(1, size(Y,2));
        
        % Fit a separate model independently for each column in Y.
        parfor i=1:size(Y,2)
            % Invoke fitlmematrix() to solve the model for this column.
            model = fitlmematrix(X, Y(:,i), Z, G, reml_args{:});
            
            % Unpack and store the results.
            B(:,i) = fixedEffects(model);
            covB(:,:,i) = model.CoefficientCovariance;
            u(:,i) = randomEffects(model);
            [psi, mse_i] = covarianceParameters(model);
            mse(1,i) = mse_i;
            sigmas(:,i) = cell2mat(cellfun(@(x)diag(x), psi, 'UniformOutput', false));
            
            if mod(i, 1000) == 0
                fprintf(1, '%0.2f%%\n', i/size(Y,2)*100);
            end
        end
        clear mse_i
    end
    
    % Compose the results.
    model = struct( ...
        'B', B, ...
        'covB', covB, ...
        'u', u, ...
        'sigmas', sigmas, ...
        'mse', mse ...
    );
elseif strcmp(method, 'swe')
    % Autodetect which columns in X correspond are also columns in Z.
    % Generates a cell array of X times Z matrices, one for each cell entry
    % in Z.  The matrix rows are columns in X.  The matrix columns are
    % columns in Z.  A matrix element of 1 (true) means the column in X is
    % also a column in Z.  The matrix is zero everywhere else.
    %
    % TODO current version of scand() does not need this functionality.
    % Consider removing completely in future version.
    %{
    XZ_maps = cell(1, length(Z));
    for i=1:length(XZ_maps)
        r = corr(X,Z{i});
        r(diag(isnan(diag(r)))) = 1;
        r(isnan(r)) = 0;
        XZ_maps{i} = r > 0.99;
    end
    clear r
    %}
    
    % Augment the matrix X with a matrix of additional fixed effect
    % covariates to account for the random effects.
    Xaug = random_to_fixed(Z, G);
    
    % Make sure we have enough degrees of freedom to account for all the
    % random effects.
    assert((size(X,2) + size(Xaug,2)) < size(X,1));
    
    % Solve the model.
    B = qr_pinv([X, Xaug]) * Y;
    [B, Baug] = deal(B(1:size(X,2), :), B(size(X,2)+1:end, :));
    
    % Approximate the random effect parameters by recombining their
    % augmented columns from B.
    u = fixed_to_random(Z, G, Baug);
    
    % Find the standard error of the fixed effects using the sandwich
    % estimator.  We need the pseudoinverse and residuals using the fixed
    % effects only.
    covB = swe_block(qr_pinv(X), Y-X*B, G, SE_only);
    
    % Compose the results.
    model = struct( ...
        'B', B, ...
        'u', u ...
    );
    if SE_only
        model.SE = covB;
    else
        model.covB = covB;
    end
elseif strcmp(method, 'mom')
    error('Not implemented.');
else
    error(['Invalid method: ' method]);
end

end

% Unpack struct (i.e. unmatched arguments) into variadic arguments.
function c = struct2opt(s)
    fname = fieldnames(s);
    fval = struct2cell(s);
    c = [fname, fval]';
    c = c(:);
end
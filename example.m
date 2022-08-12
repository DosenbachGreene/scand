% Generate example/test data for a mixed effects model.
%
% Premise: A scientist wants to test if SAT score is correlated with high
% school GPA.  She has collected data from 2300 students attending 3 high
% schools in St. Louis: Jennings (a public high school plagued by poverty),
% Gateway Science Academy (a charter school), and Burroughs (an elite
% college prep school).
%
% GPA ranges from 1.0 to 4.0 with an average of about 2.5.
% SAT scores range from 400 (lowest possible score) to 1600 (perfect).
% Both are (supposed) to be normally distributed.

%% Generate the Data

% Number of students in each school, or "cluster."
% Try playing with having an equal (balanced) vs unequal (unbalanced)
% number of students in each school.
n_jennings = 1000;
n_gateway = 800;
n_burroughs = 500;
n = n_jennings + n_gateway + n_burroughs;

% Generate GPAs.  Trim GPAs outside the range of possible values.
gpa = normrnd(2.5,0.5,n,1);
gpa(gpa < 1) = 1;
gpa(gpa > 4) = 4;

% Encode the three high schools in an n x 3 matrix.  Each high school
% occupies one column.  Values of 1 indicate a student attends the high
% school for the corresponding column.  The matrix is zero everywhere else.
school = zeros(n,3);
school(1:n_jennings, 1) = 1; % Jennings
school((n_jennings+1):(n_jennings+n_gateway), 2) = 1; % Gateway
school((n_jennings+n_gateway+1):n, 3) = 1; % Burroughs

% Also create a categorical variable to indicate school membership.
% 1 = Jennings
% 2 = Gateway
% 3 = Burroughs
G = categorical(sum(school .* [1,2,3], 2));

% Simulate some overall ground-truth relationship between SAT score and
% GPA.
intercept = 300;
slope = 200;
sat = intercept + gpa .* slope;

% Mix in the effect of school.
% Comment/uncomment to simulate random intercept, random slope, or random
% intercept + slope models.
school_intercept = [-75; 75; 300];
%school_intercept = [0; 0; 0];
school_slope = [-50; -30; 80];
%school_slope = [0; 0; 0];
sat = sat + school * school_intercept + school.*gpa * school_slope;

% Adjust the "true" slope and intercept to reflect the average cluster
% contributions.
intercept = intercept + mean(school_intercept);
slope = slope + mean(school_slope);

% The GPA of students at Burroughs is inflated.
% This makes the fixed effects only estimate biased.
% Try commenting this out and see what happens to the fixed effects only
% estimate.
gpa(school(:,3) == 1) = gpa(school(:,3) == 1) + 0.5;
gpa(gpa > 4) = 4;

% Add in homoskedastic error.
sat = sat + normrnd(0, 50, n, 1);

% Clip data to range of possible SAT scores.
sat(sat < 400) = 400;
sat(sat > 1600) = 1600;

% Store data in a closure for easy plotting of results.
plot_clusters = make_plot_fn(sat, gpa, school);

%% Fixed Effects Only Model (Without Clusters)

model_fixed = struct;

% The observations are SAT score.
Y = sat;
model_fixed.Y = Y;

% Construct fixed effects design matrix.
X = [ones(n,1), gpa];
model_fixed.X = X;

% Fit the fixed effects model.
Xpinv = pinv(X);
model_fixed.B = pinv(X)*Y;

% Scatterplot
plot_clusters(model_fixed.B);
title('Fixed Effects Without Clusters');

% Find the standard errors.
resid = model_fixed.Y - model_fixed.X * model_fixed.B;
model_fixed.SE = sqrt(diag(Xpinv*Xpinv').*sum(resid.*resid)/size(Y,1));

% Compute Standard Error for Fixed Effects Model Using Sandwich Estimator
model_fixed.SE_swe = swe_block(Xpinv, resid, G, true);
clear resid Xpinv

%% Mixed Effects Model Using REML
% Random intercept and slope.

% We could also have just called Matlab's fitlmematrix directly. Note we
% are passing in X a second time for the random effects design matrix Z.
% We reuse the intercept column for the random intercept and the slope
% column for the random slope.
model_reml = scand(X,Y,X,G, 'method', 'reml');

% Scatterplot
plot_clusters(model_reml.B, model_reml.u);
title('REML');

% Standard errors of fixed effects estimates.
model_reml.SE = covB_to_SE(model_reml.covB);

%% Fixed Effects with Clusters (Marginal Model)
% Random interecept and slope.

model_marginal = scand(X,Y,X,G, 'method', 'swe');

% Scatterplot
plot_clusters(model_marginal.B, model_marginal.u);
title('Fixed Effects With Clusters');

% Standard error with sandwich estimator.
model_marginal.SE_swe = covB_to_SE(model_marginal.covB);

% Standard error without sandwich estimator.
X = [ones(n,1), gpa, school(:,1).*1.5-0.5, school(:,2).*1.5-0.5, gpa.*(school(:,1).*1.5-0.5), gpa.*(school(:,2).*1.5-0.5)];
B = pinv(X)*Y;
model_marginal.SE = sqrt(diag(pinv(X)*pinv(X)'*sum((Y-X*B).^2)/size(Y,1)));
clear X Y B
%% Visualize comparison table.

% Observe that the fixed effects only (without clusters) estimate of the
% intercept and slope is probably biased (assuming grade inflation is not
% commented out at the top of this file).  This will happen anytime there
% is a correlation between fixed effects and random effects design matrices
% (in this case, gpa is correlated with which school a student attends.)
%
% Observe that, without using the sandwich estimator, the fixed effects
% estimates of the standard error are much too narrow, which would lead to
% inflated T-values and false positive results!  On the other hand, the SwE
% and MoM SE tend to be bigger than the REML SE, leading to low efficiency
% and false negative results.

% Fixed effects without clusters.
B = model_fixed.B(1);
SE = model_fixed.SE(1);
SE_swe = model_fixed.SE_swe(1);
boxdata_int = [B-SE_swe, B-SE, B, B+SE, B+SE_swe]; % intercept
B = model_fixed.B(2);
SE = model_fixed.SE(2);
SE_swe = model_fixed.SE_swe(2);
boxdata_slope = [B-SE_swe, B-SE, B, B+SE, B+SE_swe]; % slope
clear B SE SE_swe

% REML
B = model_reml.B(1);
SE = model_reml.SE(1);
boxdata_int = [boxdata_int; [B-SE, B-SE, B, B+SE, B+SE]];
B = model_reml.B(2);
SE = model_reml.SE(2);
boxdata_slope = [boxdata_slope; [B-SE, B-SE, B, B+SE, B+SE]];
clear B SE

% Marginal (fixed effects with clusters).
B = model_marginal.B(1);
SE = model_marginal.SE(1);
SE_swe = model_marginal.SE_swe(1);
boxdata_int = [boxdata_int; [B-SE_swe, B-SE, B, B+SE, B+SE_swe]];
B = model_marginal.B(2);
SE = model_marginal.SE(2);
SE_swe = model_marginal.SE_swe(2);
boxdata_slope = [boxdata_slope; [B-SE_swe, B-SE, B, B+SE, B+SE_swe]];
clear B SE SE_swe

% Prepare boxplot data for plotting.
boxdata_int = boxdata_int(:, [1 2 2 3 4 4 5])';
boxdata_slope = boxdata_slope(:, [1 2 2 3 4 4 5])';

% Make the plot.
f = figure; f.Position(3) = f.Position(3) * 2;
subplot(1,2,1);
boxplot(boxdata_int, 'Whisker', Inf);
xlabels = {'Fixed Only', 'REML', 'Marginal'};
xticklabels(xlabels);
ylabel('SAT Score');
title('Interecept');
clear boxdata_int
subplot(1,2,2);
h_box = boxplot(boxdata_slope, 'Whisker', Inf);
h_box = handle(h_box);
xticklabels(xlabels);
ylabel('Change in SAT Score per GPA Point');
title('Slope');
h_swe = h_box(2);
h_se = h_box(5);
legend([h_se, h_swe], 'Standard Error', 'SwE Error', 'Location', 'Northeast');
clear h_box h_se h_swe  boxdata_slope xlabels f
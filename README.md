# SCAND

Speedy Cluster Analysis for Neuroimaging Data is a reference implementation of methods for approximate but fast mixed effects linear models in Matlab. Such methods are useful in large neuroimaging data where unbiased methods such as restricted maximum likelihood (REML) may be computationally infeasible. Scand is written to be reasonably optimized, but foremost to be a readable implementation of the underlying mathematical concepts suitable as a reference for creation of other implementations. As such, scand itself may not have the features and performance to satisfy the demands of applied neuroimaging research problems.

## Work in Progress

This repository remains a work in progress. Some of the key features have not been implemented yet.

## Quick Start

Clone this git repository and open the directory in Matlab.  The function [`scand()`](./scand.m) provides a high-level interface to the different methods, which is deliberately similar to Matlab's [`fitlmematrix()`](https://www.mathworks.com/help/stats/fitlmematrix.html). Read through and try out the code in [`example.m`](./example.m).

```
git clone https://github.com/DosenbachGreene/scand.git
matlab
>> cd scand
>> help scand
 scand Approximately fit a linear mixed effects model to the data.
 
 scand (Speedy Cluster Analysis for Neuroimaging Data) is a high-level
 function for approximating a linear mixed effects model...
>> example
```

## Sandwich Estimator/Marginal Model

The sandwich estimator (SwE) method for approximate mixed effects model estimation estimates a (linear) fixed effects marginal model using ordinary least squares (OLS) and then adjusts the standard errors of the fixed effects estimates using the Huber-White sandwich estimator. The "random" cluster effects are modeled as additional fixed effects in the OLS model to account for correlations between the fixed and random effects.

The sandwich estimator method is faster than the method of moments (MoM). One tradeoff is loss of degrees of freedom when including cluster-specific terms in the OLS model. This is particularly problematic when the model contains many small clusteres, e.g. repeated measures from the same individual in a longitudinal study. The method of moments may be better suited to such problems. Another tradeoff is inability to estimate the between-cluster variance in the marginal model as can be done in method of moments.

For more details please see Bryan Guillaume's thesis work and check out the Oxford Big Data Institute's implementation of sandwich estimators at <https://www.nisox.org/Software/SwE/>.

> Guillaume B, Hua X, Thompson PM, Waldorp L, Nichols TE; Alzheimer's Disease Neuroimaging Initiative. Fast and accurate modelling of longitudinal and repeated measures neuroimaging data. Neuroimage. 2014 Jul 1;94:287-302. [doi: 10.1016/j.neuroimage.2014.03.029](https://doi.org/10.1016/j.neuroimage.2014.03.029).

Nota bene that the OLS model may or may not be augmented by additional "fixed effects" covariates to model the random cluster effects. (Scand _does_ augment by default.) If the OLS model does not contain terms to model the cluster effects then the fixed effects estimates will be biased when the fixed and cluster effects are not independent. See David Freedman's commentary, [On The So-Called "Huber Sandwich Estimator" and "Robust Standard Errors"](https://doi.org/10.1198/000313006X152207) for more about this possible source of bias.

## Method of Moments

The method of moments (MoM) uses Haeseman-Elston regression to estimate between-cluster variance from an empirical estimate of the variance covariance matrix $V$. This empirical estimate is taken as some variation of $V = \epsilon \epsilon'$ where $\epsilon$ are the residuals of the fixed effects only model. Once the between-cluster variance is estimated, it is used to generate a constrained estimate of the variance-covariance matrix $\hat{V}$ which is then used to find the fixed effect coefficients $\beta$ and random effects $u$ using the usual formulae for mixed effects models.

The matrix product $\epsilon \epsilon'$ is the second order moment of the residuals, giving rise to the name "method of moments." The methods of moment estimator is biased, but theoretically much faster to compute than the restricted maximum likelihood (REML) estimate of the between-cluster variance.

Method of moments is slower than the sandwich estimator/marginal model (SwE) approach. Furthermore, Haeseman-Elston regression can sometimes yield problematic, negative estimates for variance, which  does not occur with SwE. (These negative estimates are most likely to arise when the model is misspecified, e.g. including a random effect for slope when one does not exist.) As noted above, MoM is most useful when SwE would waste too many degrees of freedom, or when the value of the between-cluster variance (which is accounted for but not estimated by SwE) is of interest.
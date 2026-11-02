# VSOF
Code to fit variational scalar-on-function (VSOF) regression models and calculate a global test of functional effects. The file vsof.R is in the Code folder and contains code for implementing both variational scalar-on-function regression and the ZLS test. The primary functions are `vsof()` and `zls_test()`, the usage of which we describe below. Additional functions like `coef()` and `plot()` are defined for objects of class `vsof`. 

## Functions

### VSOF Model
```
vsof(fixed, se = NULL, fe = NULL, Ts = NULL, data, 
                  se.spline = NULL, fe.spline = NULL, family = gaussian(link = "identity"),
                  priors = list(Ae = 0.01, Be = 0.01, Af = c(0.01, 0.01),
                                Bf = c(0.01, 0.01), sigB = 1000),
                  tol = 1e-5, maxIter = 500, up = 100, dots = up/10,
                  verbose = TRUE)
```
`fixed` is an R `formula` object specifying the outcome and fixed effects, minimally set to `y ~ 1`. Scalar covariates can be added in the usual fashion: `y ~ x1 + x2 + ...`. To identify functional effects, use the argument `fe` which takes a one sided `formula` object of the form `~ W` where `W` is a $n\times T$ matrix of functional effects. To specify multiple functional effects, use `fe = ~ W1 + W2 + ...`. `Ts` specifies the domain of the functional effect. By default, it is taken to be the number of columns of `W` (or `W1`, `W2`, etc., in the multiple function case. `data` can be used to provide an object of class `data.frame` containing the outcome and fixed covariates.

To change the default spline specifications for the functional effect, set `fe.spline = list(df, degree, xi)` where `df` is the number of knots, `degree` is the degree of the B-spline, and `xi` is the tuning parameter on the penalty matrix. The default is `fe.spline = list(df = 8, degree = 3, xi = 0.01)`. `priors` specifies the default fixed hyperparameters. The `tol` can be altered as desired to determine convergence. `maxIter` specifies how many iterations the model is allowed to run before terminating should convergence not be reached. `verbose = TRUE` prints updates of the model's progress, with updates depending on `up` and `dots`.

#### Experimental Features
The argument `se` can be used to specify smoothed, non-functional effects. The argument `family` can take both `gaussian(link = "identity")` (as shown) and `binomial(link = "probit")` for binary outcomes. Other link functions are not currently supported. The smoothed effect and binary outcome are experimental and not fully vetted. Meyer and Cheng (2026) only considers Gaussian outcomes and functional effects.

#### Output
`vsof()` returns an object of class `vsof`. Valid class functions include `plot()`, `coef()`, and `zls_test()` (see below). The object itself contains the variational estimates of the mean in `out$VBEstimates`, the transformed functional effects in `out$funcEffect`, variance estimates in `out$varComps`, among other details related to the model fit. Convergence details can be found in `out$algOut` including the number of iterations (`iter`), the ELBO (`logpi`), and changes in the ELBO (`lDelta`). The vector required by Theorem 1 of manuscript is `out$Ctf`.

### ZLS Test
```
zls_test(model, test = NULL)
```
`model` is of class `vsof`. It performs the ZLS test for variational functional effects and returns an object of class `htest`. When a model has multiple functional effects, all tests will be returned if `test = NULL`. The argument `test` can be used to extract a specific test when multiple functional effects are fit. The order the tests appear in is the order implied by the one-sided `formula` supplied to the `fe` argument. For example, if the `fe` argument has the form `fe = ~ W1 + W2 + W3`, then `zls_test(model, test = 1)` extracts the test for `W1`, `zls_test(model, test = 2)` extracts the test for `W2`, and `zls_test(model, test = 3)` extracts the test for `W3`.

## Simulations folder
Contains sample code for the parallelized run of simulation settings described by Meyer and Cheng (2026). The file bayesian_helpers.R in the Code folder is required to implement the fully Bayesian models in the simulation file.

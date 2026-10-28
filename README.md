# VSOF
Code to fit variational scalar-on-function (VSOF) regression models and calculate a global test of functional effects. The file vsof.R is in the Code folder and contains code for implementing both variational scalar-on-function regression and the ZLS test. The primary functions are `vsof()` and `zls_test()`, the usage of which we describe below. Additional functions like `coef()` and `plot()` are defined for objects of class `vsof`. 

## Functions

### VSOF model
```
vsof(fixed, se = NULL, fe = NULL, Ts = NULL, data, 
                  se.spline = NULL, fe.spline = NULL, family = gaussian(link = "identity"),
                  priors = list(Ae = 0.01, Be = 0.01, Af = c(0.01, 0.01), # edit defaults
                                Bf = c(0.01, 0.01), sigB = 1000),
                  tol = 1e-5, maxIter = 500, up = 100, dots = up/10,
                  verbose = TRUE)
```
`fixed` is an R `formula` object specifying the outcome and fixed effects, minimally set to `y ~ 1`. To identify functional effects, use the argument `fe` which takes a one sided `formula` object of the form `~ W` where `W` is a $n\times T$ matrix of functional effects. To specify multiple functional effects, use `fe = ~ W1 + W2`, etc.

### ZLS test
```
zls_test(model, test = NULL)
```
`model` is of class `vsof`. It performs the ZLS test for variational functional effects and returns an object of class `htest`. When a model has multiple functional effects, all tests will be returned if `test = NULL`. The argument `test` can be used to extract a specific test when multiple functional effects are fit. The order the effects appear in the `fe` one-sided `formula` determine the numerical order of the test, 1, 2, etc.

## Simulations folder
Contains sample code for the parallelized run of simulation settings described by Meyer and Cheng (2026). The file bayesian_helpers.R in the Code folder is required to implement the fully Bayesian models in the simulation file.

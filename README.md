# VSOF
Code to fit variational scalar-on-function (VSOF) regression models and calculate a global test of functional effects.

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

### ZLS test
```
zls_test(...)
```

## Simulations folder
Contains sample code for the parallelized run of simulation settings described by Meyer and Cheng (2026).

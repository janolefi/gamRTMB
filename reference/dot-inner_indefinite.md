# Diagnose a non-finite marginal objective as non-concavity

A non-finite marginal objective is usually read as the response leaving
the family's support, and the error message used to say so. That is the
wrong diagnosis for a whole class of models, and a misleading one: the
Laplace approximation needs the log determinant of the inner Hessian, so
an *indefinite* Hessian produces exactly the same symptom with the data
entirely inside the support.

## Usage

``` r
.inner_indefinite(obj, design, method, famname)
```

## Arguments

- obj:

  The `MakeADFun` object, evaluated at its starting values.

- design:

  The design object.

- method:

  `"REML"` or `"ML"`.

- famname:

  The family's name, for the message.

## Value

A sentence describing the negative curvature, or `NULL` if the Hessian
is unavailable or positive definite.

## Details

The Box-Cox power exponential is the case in hand. Its log-likelihood is
concave in the four intercepts alone – an intercept-only fit converges
in half a second – but adding a single covariate column to `mu` puts
three negative eigenvalues into the inner Hessian, every one of them a
direction mixing that column with the `sigma`, `nu` and `tau`
intercepts. No starting value repairs it: sweeping `tau` from 2 to 9 and
`nu` from 1 to 2.5 never gets below two negative directions. It is a
property of the family's parameterisation, not of the start.

Under `"REML"` the negative curvature lands in `beta`, which is declared
random and so passes through the inner solve carrying no prior to
convexify it. The penalized blocks are not the problem – their Gaussian
prior leaves them comfortably positive definite. Hence the suggestion of
a basis with no null space, which is what puts those columns under a
penalty.

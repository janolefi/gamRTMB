# Per-distribution corrections to the name-based defaults

`links` overrides the dictionary, `fixed` marks arguments that are known
data rather than parameters, `modelled` pins the subset.

## Usage

``` r
.family_overrides
```

## Format

An object of class `list` of length 24.

## Details

The ambiguous names, and why each must be resolved per distribution:

- alpha:

  a positive shape in frechet/llogis/kumar, an unconstrained skewness in
  skewnorm/skewnorm2/sn.

- nu:

  an unconstrained Box-Cox power in bccg/bcpe/bct/gengamma, a positive
  power in powerexp, positive in combinom.

- theta:

  a positive rate in bell, a direction vector in vmf2.

- size:

  the overdispersion *parameter* in nbinom2-type densities, the known
  number of trials in binomial-type ones. The most consequential
  distinction in the table.

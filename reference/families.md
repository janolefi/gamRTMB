# Available families

Every density
[`fam()`](https://janolefi.github.io/gamRTMB/reference/fam.md) can turn
into a family: all of RTMBdist that is a univariate regression family,
plus the standard RTMB densities.

## Usage

``` r
families(pattern = NULL)
```

## Arguments

- pattern:

  Optional regular expression to filter family names.

## Value

A data frame with one row per family: its name, the modelled parameters
with their links, any fixed arguments that must be supplied from the
data, the response support, whether quantile residuals are available
(i.e. whether a CDF exists), and which package the density comes from.

## See also

[`fam()`](https://janolefi.github.io/gamRTMB/reference/fam.md)

## Examples

``` r
head(families(), 10)
#>       family                              parameters needs    support residuals
#> 1       bccg          mu/log, sigma/log, nu/identity       continuous      TRUE
#> 2       bcpe mu/log, sigma/log, nu/identity, tau/log       continuous      TRUE
#> 3        bct mu/log, sigma/log, nu/identity, tau/log       continuous      TRUE
#> 4       bell                               theta/log          lattice      TRUE
#> 5      bell2                                  mu/log          lattice      TRUE
#> 6       beta                  shape1/log, shape2/log       continuous      TRUE
#> 7      beta2                       mu/logit, phi/log       continuous      TRUE
#> 8  betabinom                  shape1/log, shape2/log  size    lattice     FALSE
#> 9  betaprime                  shape1/log, shape2/log       continuous      TRUE
#> 10     binom                              prob/logit  size    lattice      TRUE
#>    quantiles   source
#> 1       TRUE RTMBdist
#> 2       TRUE RTMBdist
#> 3       TRUE RTMBdist
#> 4       TRUE RTMBdist
#> 5       TRUE RTMBdist
#> 6       TRUE RTMBdist
#> 7       TRUE RTMBdist
#> 8      FALSE RTMBdist
#> 9       TRUE RTMBdist
#> 10      TRUE     RTMB
families("beta")
#>       family                                            parameters needs
#> 1       beta                                shape1/log, shape2/log      
#> 2      beta2                                     mu/logit, phi/log      
#> 3  betabinom                                shape1/log, shape2/log  size
#> 4  betaprime                                shape1/log, shape2/log      
#> 5     oibeta                 shape1/log, shape2/log, oneprob/logit      
#> 6    oibeta2                      mu/logit, phi/log, oneprob/logit      
#> 7     zibeta                shape1/log, shape2/log, zeroprob/logit      
#> 8    zibeta2                     mu/logit, phi/log, zeroprob/logit      
#> 9    zoibeta shape1/log, shape2/log, zeroprob/logit, oneprob/logit      
#> 10  zoibeta2      mu/logit, phi/log, zeroprob/logit, oneprob/logit      
#>       support residuals quantiles   source
#> 1  continuous      TRUE      TRUE RTMBdist
#> 2  continuous      TRUE      TRUE RTMBdist
#> 3     lattice     FALSE     FALSE RTMBdist
#> 4  continuous      TRUE      TRUE RTMBdist
#> 5       mixed      TRUE     FALSE RTMBdist
#> 6       mixed      TRUE     FALSE RTMBdist
#> 7       mixed      TRUE     FALSE RTMBdist
#> 8       mixed      TRUE     FALSE RTMBdist
#> 9       mixed      TRUE     FALSE RTMBdist
#> 10      mixed      TRUE     FALSE RTMBdist
```

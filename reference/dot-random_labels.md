# Label every coefficient the inner problem solves over

In the order
[`RTMB::MakeADFun()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html)
lays them out, which is the order of the parameter list: all of `beta`,
then all of `b`. Under `"ML"` only `b` is random, so only those are
named.

## Usage

``` r
.random_labels(design, method)
```

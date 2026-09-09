# gamRTMB: distributional regression with mgcv smooths and RTMB

Fits GAMLSS-style models in which every parameter of a distribution may
carry smooth terms. The pieces are borrowed rather than rebuilt: mgcv
constructs the bases and penalties,
[`mgcv::smooth2random()`](https://rdrr.io/pkg/mgcv/man/smooth2random.html)
turns the penalized coefficients into iid Gaussian random effects, RTMB
supplies automatic differentiation and the Laplace approximation, and
RTMBdist supplies the log-densities in each distribution's own native
parameterisation.

## Entry points

[`gamRTMB()`](https://janolefi.github.io/gamRTMB/reference/gamRTMB.md)
fits a model,
[`fam()`](https://janolefi.github.io/gamRTMB/reference/fam.md) builds a
family object from any suitable density,
[`families()`](https://janolefi.github.io/gamRTMB/reference/families.md)
lists what is available, and
[`edf()`](https://janolefi.github.io/gamRTMB/reference/edf.md) reports
effective degrees of freedom per smooth.

## References

Wood, S. N. (2011) Fast stable restricted maximum likelihood and
marginal likelihood estimation of semiparametric generalized linear
models. *JRSS-B* 73(1), 3–36.

Wood, S. N. and Fasiolo, M. (2017) A generalized Fellner-Schall method
for smoothing parameter optimization with application to Tweedie
location, scale and shape models. *Biometrics* 73(4), 1071–1081.

## See also

Useful links:

- <https://github.com/janolefi/gamRTMB>

- <https://janolefi.github.io/gamRTMB/>

- Report bugs at <https://github.com/janolefi/gamRTMB/issues>

## Author

**Maintainer**: Jan-Ole Fischer <jan-ole.fischer@mailbox.org>
([ORCID](https://orcid.org/0009-0004-1556-9053))

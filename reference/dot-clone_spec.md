# Copy a basis specification onto another term

Local stand-in for the unexported `mgcv:::clone.smooth.spec`: the first
member of an `id` group defines the basis, and each later term gets that
definition back with its own identity (term names, label, `by`, `xt`).

## Usage

``` r
.clone_spec(base, spec)
```

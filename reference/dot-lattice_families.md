# Lattice (integer-supported) families

Whether a response is discrete cannot be read off a density: the
argument names say nothing about it, and probing the CDF numerically is
unreliable because implementations differ in whether they floor a
non-integer argument. So it is declared, as in LaMa's `pseudo_res()`,
and [`fam()`](https://janolefi.github.io/gamRTMB/reference/fam.md) takes
a `support` override for anything not listed.

## Usage

``` r
.lattice_families
```

## Format

An object of class `character` of length 29.

## Details

Every family here is supported on the non-negative integers, which is
what lets
[`residuals.gamRTMB()`](https://janolefi.github.io/gamRTMB/reference/residuals.gamRTMB.md)
take `F(y - 1) = 0` at `y = 0` rather than evaluating the CDF below its
support.

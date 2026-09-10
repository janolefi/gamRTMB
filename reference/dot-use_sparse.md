# Is this smooth's penalty worth keeping sparse?

Two ways to qualify. A smooth with a single penalty qualifies if that
penalty is big and sparse enough that `smooth2random`'s rotation would
actually cost something – a 10-coefficient P-spline penalty is
technically banded, but nothing is gained by treating it specially and
the well-travelled route is the safer one.

## Usage

``` r
.use_sparse(sm, mode)
```

## Arguments

- sm:

  A `smoothCon` object.

- mode:

  `"auto"`, `"always"` or `"never"`.

## Details

A smooth carrying an `L` matrix qualifies outright, whatever its size.
`L` is mgcv's way of saying that several penalty matrices combine into
one through fewer smoothing parameters than there are matrices, as
\\\lambda = \exp(L\theta)\\. `smooth2random` cannot represent that at
all – it needs one variance per penalized block – so such a smooth has
nowhere else to go. The SPDE smooth is the case in hand: three finite
element matrices, two parameters.

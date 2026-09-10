# A length scale for the meshed domain, for starting values

The larger side of the bounding box of the mesh nodes that actually
carry basis weight, which is the region the data occupy rather than the
outer extension. A one-dimensional mesh stores its knots as a plain
vector, and its degree-2 basis has one more function than it has knots,
so both the shape and the length of `loc` have to be taken as they come.

## Usage

``` r
.spde_extent(mesh, X)
```

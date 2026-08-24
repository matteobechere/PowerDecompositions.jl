# PowerDecompositions.jl

A Julia package for numerically finding and verifying decompositions of polynomials as sums of powers of forms (homogeneous polynomials).

## Installation

To install the package, in Julia, enter the package manager by pressing `]` and then type

`add https://github.com/matteobechere/PowerDecompositions.jl`.

## Basic usage

```
using PowerDecompositions
using DynamicPolynomials

@polyvar x y

F = x^3 + y^3
result = decompose_powers(F, 3)

print_decomposition(result)
```

## Output

The function `decompose_powers` returns a named tuple called `result`, describing a decomposition

$$
F \approx \sum_{i=1}^{r} \omega_i q_i^d,
$$

where each $q_i$ is a form of degree $h=\frac{\deg(F)}{d}$.

In particular, the $q_i$ are linear forms only when $\deg(F)=d$.

The `result` contains:

- `verified`: a Boolean indicating whether the decomposition was successfully verified;
- `field`: the field used for the decomposition, either `:real` or `:complex`;
- `power`: the exponent $d$;
- `rank`: the number $r$ of summands returned;
- `weights`: a vector containing the scalar coefficients $\omega_1,\ldots,\omega_r$. When `absorb_weights=true`, scaling is absorbed into the forms whenever possible;
- `forms`: a vector containing the degree-$h$ forms $q_1,\ldots,q_r$;
- `residual`: the relative coefficient norm of the difference between $F$ and the reconstructed polynomial. A small `residual` indicates an accurate decomposition;
- `message`: a human-readable description of the result.

If no verified decomposition is found, `verified` is `false`, `rank` is `nothing`, `weights` and `forms` are empty, and `residual` is `Inf`.

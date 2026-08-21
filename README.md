# PowerDecompositions.jl

A Julia package for numerically finding and verifying decompositions of polynomials as sums of powers of linear forms.

## Installation

To install the package, in Julia, enter the package manager by pressing `]` and then type

`add https://github.com/matteobechere/PowerDecompositions.jl`.

## Basic usage

using PowerDecompositions
using DynamicPolynomials

@polyvar x y

F = x^3 + y^3
result = decompose_powers(F, 3)

print_decomposition(result)

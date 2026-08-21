# PowerDecompositions.jl

A Julia package for numerically finding and verifying decompositions of polynomials as sums of powers of linear forms.

## Installation

From the Julia REPL:

```julia
using Pkg
Pkg.add(url="https://github.com/matteobechere/PowerDecompositions")

## Basic usage

using PowerDecompositions
using DynamicPolynomials

@polyvar x y

F = x^3 + y^3
result = decompose_powers(F, 3)

print_decomposition(result)
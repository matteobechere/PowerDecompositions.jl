# Return the degree of a nonzero homogeneous polynomial
function homogeneous_degree(F)
    ts = terms(F)
    isempty(ts) && throw(ArgumentError("F must be nonzero"))
    ds = unique(degree(monomial(t)) for t in ts)
    length(ds) == 1 || throw(ArgumentError("F must be homogeneous"))
    only(ds)
end

# Construct the affine family of degree-d preimages of F
function preimage_family(F, d::Int; X=collect(variables(F)))
    d > 0 || throw(ArgumentError("d must be positive"))
    X = collect(X)
    isempty(X) && throw(ArgumentError("at least one variable is required"))
    all(x -> x isa MultivariatePolynomials.AbstractVariable, X) ||
        throw(ArgumentError("every element of X must be a polynomial variable"))
    allunique(X) || throw(ArgumentError("X must not contain duplicate variables"))
    all(x -> x in X, variables(F)) ||
        throw(ArgumentError("X must contain every variable of F"))

    D = homogeneous_degree(F)
    D % d == 0 || throw(ArgumentError("degree(F) must be divisible by d"))
    h = D ÷ d
    Xbeta = monomials(X, h)

    @polyvar Y[1:length(Xbeta)]
    Ymons = monomials(Y, d)
    images = [
        prod(Xbeta[i]^a for (i, a) in enumerate(exponents(m)))
        for m in Ymons
    ]

    fibres = Dict{eltype(images),Vector{Int}}()
    for (i, image) in enumerate(images)
        push!(get!(fibres, image, Int[]), i)
    end

    G = sum(
        coefficient(t) * Ymons[first(fibres[monomial(t)])]
        for t in terms(F);
        init=zero(polynomial(first(Ymons))),
    )
    kernel = [
        Ymons[j] - Ymons[first(I)]
        for I in values(fibres) for j in Iterators.drop(I, 1)
    ]

    @polyvar tau[1:length(kernel)]
    P = G + sum(
        tau[i] * kernel[i] for i in eachindex(kernel);
        init=zero(G),
    )
    (; P, G, kernel, tau, Y, Xbeta, h)
end
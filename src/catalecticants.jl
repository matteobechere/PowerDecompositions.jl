# Treat non-Y variables as coefficients and index them by Y exponents
function y_coefficients(P, Y)
    out = Dict{Tuple{Vararg{Int}},Any}()
    for t in terms(P)
        m = monomial(t)
        alpha = Tuple(degree(m, y) for y in Y)
        c = coefficient(t)
        for x in variables(m)
            x in Y || (c *= x^degree(m, x))
        end
        out[alpha] = get(out, alpha, 0) + c
    end
    out
end

multinomial(alpha) =
    factorial(big(sum(alpha))) ÷ prod(factorial(big(a)) for a in alpha)

# Construct the catalecticant of a homogeneous polynomial
function catalecticant(P, Y; split=nothing)
    ydegrees = unique(
        sum(degree(monomial(t), y) for y in Y) for t in terms(P)
    )
    length(ydegrees) == 1 ||
        throw(ArgumentError("the preimage must be homogeneous in Y"))
    d = only(ydegrees)
    a = isnothing(split) ? d ÷ 2 : Int(split)
    0 <= a <= d || throw(ArgumentError("split must be between 0 and $d"))

    rows, cols = monomials(Y, a), monomials(Y, d - a)
    coeffs = y_coefficients(P, Y)
    prototype = (1 // 1) * polynomial(P)
    C = Matrix{typeof(prototype)}(undef, length(rows), length(cols))
    for i in eachindex(rows), j in eachindex(cols)
        gamma = exponents(rows[i]) .+ exponents(cols[j])
        C[i, j] = get(coeffs, Tuple(gamma), 0) / multinomial(gamma)
    end
    (; C, rows, cols)
end

constant_value(p::Number) = p
function constant_value(p)
    cs = coefficients(p)
    isempty(cs) ? 0 : only(cs)
end

numeric_value(z, ::Type{Float64}) = Float64(real(z))
numeric_value(z, ::Type{ComplexF64}) = ComplexF64(z)

function evaluate_matrix(C, tau, values, ::Type{T}) where {T}
    M = Matrix{T}(undef, size(C))
    for i in axes(C, 1), j in axes(C, 2)
        q = isempty(tau) ? C[i, j] : subs(C[i, j], tau => values)
        M[i, j] = numeric_value(constant_value(q), T)
    end
    M
end

# Write the catalecticant as C(tau) = C0 + sum(tau[i] * A[i])
function affine_catalecticant(C, tau, ::Type{T}) where {T}
    n = length(tau)
    z = zeros(T, n)
    C0 = evaluate_matrix(C, tau, z, T)
    A = Matrix{T}(undef, length(C0), n)
    for j in 1:n
        e = zeros(T, n)
        e[j] = one(T)
        A[:, j] = vec(evaluate_matrix(C, tau, e, T) - C0)
    end
    C0, A
end
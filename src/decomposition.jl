
# Converts a symbolic polynomial to floating-point coefficients
numeric_polynomial(P, ::Type{T}) where {T} = sum(
    numeric_value(coefficient(t), T) * monomial(t) for t in terms(P);
    init=zero(T) * first(monomials(P)),
)

# Measures the relative error of the reconstructed polynomial w.r.t. the original one
function polynomial_residual(F, approximation)
    target = ComplexF64.(coefficients(F))
    err = ComplexF64.(coefficients(approximation - F))
    norm(err) / max(norm(target), eps(Float64))
end

# Tests whether a complex number is numerical
is_approximately_real(z; atol, rtol) =
    abs(imag(z)) <= atol + rtol * abs(z)

# Calls TensorDec.decompose to perform a Waring decomposition
function tensor_decomposition_attempt(
    F,
    numeric_preimage,
    preimage,
    d,
    field,
    rank_value,
    mode;
    tol,
    atol,
    rtol,
    absorb_weights,
    verbose,
)
    answer = try
        rkf = mode == :catalecticant_rank ?
            AlgebraicSolvers.cst_rkf(rank_value) :
            AlgebraicSolvers.eps_rkf(1e-6)

        TensorDec.decompose(
            numeric_preimage,
            rkf,
            :Random;
            verbose,
        )
    catch err
        return (state=:tensor_solver_failed, message=sprint(showerror, err))
    end

    try
        weights, points = collect(answer[1]), answer[2]
        if field == :real
            if any(!is_approximately_real(z; atol, rtol) for z in weights) ||
               any(!is_approximately_real(z; atol, rtol) for z in points)
                return (
                    state=:complex_result,
                    message="the tensor solver returned a complex result",
                )
            end
            weights, points = real.(weights), real.(points)
        end

        used_Y = collect(variables(numeric_preimage))
        used_Xbeta = [
            preimage.Xbeta[findfirst(==(y), preimage.Y)] for y in used_Y
        ]
        forms = Any[
            sum(points[j, i] * used_Xbeta[j] for j in axes(points, 1))
            for i in axes(points, 2)
        ]

        if absorb_weights
            for i in eachindex(weights)
                w = weights[i]
                root = if field == :complex
                    ComplexF64(w)^(1 / d)
                elseif isodd(d)
                    sign(w) * abs(w)^(1 / d)
                elseif w >= -(atol + rtol * abs(w))
                    max(w, 0)^(1 / d)
                else
                    continue
                end
                forms[i] *= root
                weights[i] = one(w)
            end
        end

        approximation = sum(
            weights[i] * forms[i]^d for i in eachindex(forms);
            init=zero(ComplexF64) * F,
        )
        residual = polynomial_residual(F, approximation)
        residual <= tol || return (
            state=:residual_too_large,
            message="relative residual $residual exceeds tolerance $tol",
            residual=residual,
        )

        (
            state=:verified,
            message="numerically verified decomposition found",
            weights=weights,
            forms=forms,
            residual=residual,
        )
    catch err
        (state=:tensor_result_invalid, message=sprint(showerror, err))
    end
end

"""
    decompose_powers(F, d; X=variables(F), field=:auto, ...)

Numerically find and verify

    F = sum(weights[i] * forms[i]^d).

The input polynomial must have real coefficients. Set `field=:complex` to
allow complex decomposition forms for a real input polynomial.

"""
function decompose_powers(
    F,
    d::Int;
    X=collect(variables(F)),
    field::Symbol=:auto,
    split=nothing,
    attempts::Int=8,
    max_rank_iterations::Int=200,
    tol::Real=1e-9,
    atol::Real=1e-12,
    rtol::Real=1e-9,
    absorb_weights::Bool=true,
    verbose::Bool=false,
)
    all(c -> c isa Number, coefficients(F)) ||
        throw(ArgumentError("numeric coefficients are required"))
    attempts > 0 || throw(ArgumentError("attempts must be positive"))
    max_rank_iterations > 0 ||
        throw(ArgumentError("max_rank_iterations must be positive"))
    tol > 0 || throw(ArgumentError("tol must be positive"))
    field in (:auto, :real, :complex) ||
        throw(ArgumentError("field must be :auto, :real, or :complex"))

    input_is_real = all(isreal, coefficients(F))
    input_is_real || throw(ArgumentError(
        "complex input coefficients are not supported; F must have real " *
        "coefficients. Use field=:complex only to allow complex decomposition " *
        "forms for a real input polynomial.",
    ))
    chosen_field = field == :auto ? :real : field

    preimage = preimage_family(F, d; X)
    rank_result = minimum_preimage_rank(
        preimage.P,
        preimage.Y,
        preimage.tau;
        split,
        real_only=true,
        attempts,
        max_iterations=max_rank_iterations,
        tol,
    )

    if d == 2 && preimage.h == 1 && rank_result.rank > 1
        throw(ArgumentError(
            "quadratic inputs of rank greater than one are not supported",
        ))
    end

    failures = String[]
    P = rank_result.preimage
    numeric_preimage = numeric_polynomial(P, Float64)

    for mode in (:catalecticant_rank, :automatic_rank)
        for _ in 1:attempts
            outcome = tensor_decomposition_attempt(
                F,
                numeric_preimage,
                preimage,
                d,
                chosen_field,
                rank_result.rank,
                mode;
                tol,
                atol,
                rtol,
                absorb_weights,
                verbose,
            )
            if outcome.state == :verified
                return (
                    verified=true,
                    field=chosen_field,
                    power=d,
                    rank=length(outcome.forms),
                    weights=outcome.weights,
                    forms=outcome.forms,
                    residual=outcome.residual,
                    message=outcome.message,
                )
            end
            push!(failures, outcome.message)
        end
    end

    (
        verified=false,
        field=chosen_field,
        power=d,
        rank=nothing,
        weights=Any[],
        forms=Any[],
        residual=Inf,
        message=isempty(failures) ?
            "no decomposition found" : join(unique(failures), "; "),
    )
end

function print_decomposition(result)
    println("Verified: ", result.verified)
    result.verified || return println(result.message)

    println("Field: ", result.field)
    println("Rank: ", result.rank)
    println("Residual: ", result.residual)

    println("\nDecomposition:")
    for i in eachindex(result.forms)
        println(
            "  ",
            result.weights[i],
            " * (",
            result.forms[i],
            ")^",
            result.power,
        )
    end
end
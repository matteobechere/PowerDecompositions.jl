# Measures how far is a matrix M from having rank at most r
function rank_residual(M, r)
    s = svdvals(M)
    r >= length(s) && return 0.0
    Float64(norm(s[(r + 1):end]) / max(norm(s), eps(Float64)))
end

# Computes rank r approximation on M
function rank_projection(M, r)
    r == 0 && return zeros(eltype(M), size(M))
    S = svd(M)
    S.U[:, 1:r] * Diagonal(S.S[1:r]) * S.Vt[1:r, :]
end

"""
Searches for parameter values minimizing catalecticant rank by alternating projections between
the affine family C(tau) and the set of matrices of rank at most r.

For r=0,...,maxrank
    repeat max_iteration times
        Takes random element in the affine family -> Projects to matrices of rank r -> Projects this
        projection to the affine family
    r=r+1
"""
function minimum_preimage_rank(
    P,
    Y,
    tau;
    split=nothing,
    real_only=true,
    attempts=8,
    max_iterations=200,
    tol=1e-9,
)
    cat = catalecticant(P, Y; split)
    T = real_only ? Float64 : ComplexF64
    C0, A = affine_catalecticant(cat.C, tau, T)
    nparams = length(tau)
    maxrank = min(size(C0)...)
    rank_tol = tol / 100
    # First start from the chosen preimage G
    scale = nparams == 0 ? 0.0 :
        norm(C0) / max(norm(A), eps(Float64))

    for r in 0:maxrank
        best_residual = Inf
        best_values = zeros(T, nparams)

        for attempt in 1:(r == maxrank ? 1 : attempts)
            values = if attempt == 1 || nparams == 0
                zeros(T, nparams)
            elseif real_only # Then, try randomized sets of parameters
                scale .* randn(nparams)
            else
                scale .* (randn(nparams) .+ im .* randn(nparams)) ./ sqrt(2)
            end

            for _ in 0:max_iterations
                M = C0 + reshape(A * values, size(C0))
                residual = rank_residual(M, r)
                if residual < best_residual
                    best_residual = residual
                    best_values = copy(values)
                end
                residual <= rank_tol && break
                nparams == 0 && break

                low_rank_M = rank_projection(M, r)
                values = A \ vec(low_rank_M - C0)
            end

            best_residual <= rank_tol && break
        end

        if best_residual <= rank_tol
            Q = subs(P, tau => best_values)
            return (
                rank=r,
                preimage=Q,
            )
        end
    end

    error("internal error: full catalecticant rank was not reached")
end
using PowerDecompositions
using DynamicPolynomials
using Random


@polyvar x y z

"""
Run and print one example without stopping the rest of the file if the
solver does not find a decomposition.
"""
function run_example(name, F, d; seed, kwargs...)
    println("\n", "="^78)
    println(name)
    println("F = ", F)
    println("Requested power d = ", d)
    println("Degree of the returned forms = degree(F)/d")

    Random.seed!(seed)
    result = try
        decompose_powers(
            F,
            d;
            attempts=16,
            max_rank_iterations=400,
            kwargs...,
        )
    catch err
        println("The example raised an error: ", sprint(showerror, err))
        return nothing
    end

    print_decomposition(result)
    result
end

results = Dict{String,Any}()

# 1. One cube of a linear form.
results["rank_one_cubic"] = run_example(
    "1. Rank-one binary cubic",
    (x + 2y)^3,
    3;
    seed=20260822,
)

# 2. A sum of two cubes of linear forms.
results["binary_cubic"] = run_example(
    "2. Binary cubic: sum of two cubes",
    (x + y)^3 + (x - 2y)^3,
    3;
    seed=20260823,
)

# 3. Three variables and two cubic summands.
results["ternary_cubic"] = run_example(
    "3. Ternary cubic",
    (x + y + z)^3 + (2x - y + z)^3,
    3;
    seed=20260824,
)

# 4. A binary quartic written as fourth powers of linear forms.
results["binary_quartic"] = run_example(
    "4. Binary quartic",
    (x + y)^4 + 2 * (x - y)^4,
    4;
    seed=20260825,
)

# 5. A ternary quartic with two summands.
results["ternary_quartic"] = run_example(
    "5. Ternary quartic",
    (x + y + z)^4 + (x - y + z)^4,
    4;
    seed=20260826,
)

# 6. A binary sextic written as sixth powers of linear forms.
results["binary_sextic"] = run_example(
    "6. Binary sextic",
    (x + 2y)^6 + (2x - y)^6,
    6;
    seed=20260827,
)

# 7. Here degree(F)=6 and d=3, so the returned forms have degree two.
results["cubes_of_quadratics"] = run_example(
    "7. Degree-six polynomial: cubes of quadratic forms",
    (x^2 + y^2)^3 + (x * y)^3,
    3;
    seed=20260828,
)

# 8. Here degree(F)=8 and d=4, so the returned form is quadratic.
results["fourth_powers_of_quadratics"] = run_example(
    "8. Degree-eight polynomial: fourth powers of quadratic forms",
    (x^2 + x * y + y^2)^4,
    4;
    seed=20260829,
)

# 9. With absorb_weights=false, the scalar coefficients are kept separate from
#    the returned forms instead of being absorbed into them.
results["explicit_weights"] = run_example(
    "9. Cubic with explicit scalar weights",
    3 * (x + y)^3 - 2 * (x - 2y)^3,
    3;
    seed=20260830,
    absorb_weights=false,
)

# 10. A negative coefficient of an even power can be absorbed by allowing a
#     complex-valued form. The input polynomial itself still has real
#     coefficients.
results["complex_field"] = run_example(
    "10. Complex decomposition allowed",
    x^4 - y^4,
    4;
    seed=20260831,
    field=:complex,
)

println("\n", "="^78)
println("Finished ", length(results), " examples.")
println(
    "Verified decompositions: ",
    count(result -> !isnothing(result) && result.verified, values(results)),
    "/",
    length(results),
)

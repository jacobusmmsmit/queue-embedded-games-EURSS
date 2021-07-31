using Plots

include("helper_functions.jl")

begin
    λ1 = 1.0
    λ2 = 1.0
    α = 1.0
    μ(α) = 1 / (1.5 * (α))
    ν(α) = 1 / (0.7 * (α))
    qs = 0:0.01:1
    parameters = (λ1, λ2, α, μ, ν)
end

# function pstar(q, λ1, λ2, α, μ, ν)
#     p = (λ1*Complex(μ(α)*(ν(α) - λ2*q))^(1/2) - μ(α)*Complex(μ(α)*(ν(α) - λ2*q))^(1/2) - ν(α)*Complex(μ(α)*(ν(α) - λ2*q))^(1/2) - λ1*ν(α) + 2*μ(α)*ν(α) + λ2*q*Complex(μ(α)*(ν(α) - λ2*q))^(1/2) + λ1*λ2*q - 2*λ2*μ(α)*q)/(λ1*μ(α) - λ1*ν(α) + λ1*λ2*q) |> real
#     return clamp(p, 0, 1)
# end


# function best_responded_cost(p, λ1, λ2, α, μ, ν)
#     cost(p, pstar(p, λ2, λ1, α, μ, ν), parameters...)
# end

# function ∂best_responded_cost_∂p(p, λ1, λ2, α, μ, ν)
#     ForwardDiff.derivative(ϕ -> best_responded_cost(ϕ, λ1, λ2, α, μ, ν), p)
# end

# using BenchmarkTools
# brc = @benchmark best_responded_cost(0.5, λ1, λ2, α, μ, ν)
# abrc = @benchmark analytic_best_responded_cost(0.5, λ1, λ2, α, μ, ν)


# function analytic_best_responded_cost(p, λ1, λ2, α, μ, ν)
#     real(- (p - 1)/(μ(α) - λ1 + λ1*p) - p/(λ1*p - ν(α) + (λ2*(λ2*Complex(μ(α)*(ν(α) - λ1*p))^(1/2) - μ(α)*Complex(μ(α)*(ν(α) - λ1*p))^(1/2) - ν(α)*Complex(μ(α)*(ν(α) - λ1*p))^(1/2) - λ2*ν(α) + 2*μ(α)*ν(α) + λ1*p*Complex(μ(α)*(ν(α) - λ1*p))^(1/2) + λ1*λ2*p - 2*λ1*μ(α)*p))/(λ2*μ(α) - λ2*ν(α) + λ1*λ2*p)))
# end

# function diff_best_responded_cost(p, λ1, λ2, α, μ, ν)
#     squareroot_factor = real(Complex(μ(α)*(ν(α) - λ1*p))^(1/2))
#     A = (
#         (λ1*(p - 1)) / (μ(α) - λ1 + λ1*p)^2 -
#         (μ(α) - ν(α) + λ1*p)/
#         ((squareroot_factor - ν(α) + λ1*p)*(λ2 - μ(α) - ν(α) + λ1*p)) -
#         1/(μ(α) - λ1 + λ1*p)
#     )  
#     B = (
#         p * λ1 *(
#             4*real(Complex(μ(α)*(ν(α) - λ1*p))^(3/2)) -
#             μ(α)^3 + 
#             λ2*μ(α)^2 - 
#             2*λ1^2*p^2*squareroot_factor -
#             4*μ(α)^2*(ν(α) - λ1*p) +
#             2*μ(α)^2*squareroot_factor -
#             2*ν(α)^2*squareroot_factor -
#             λ1*μ(α)*p*(ν(α) - λ1*p) +
#             4*λ1*ν(α)*p*squareroot_factor +
#             λ2*μ(α)*(ν(α) - λ1*p) +
#             μ(α)*ν(α)*(ν(α) - λ1*p) -
#             2*λ2*μ(α)*squareroot_factor
#         )
#     )

#     C = (
#         2*squareroot_factor *
#         (
#             squareroot_factor -
#             ν(α) +
#             λ1*p
#         )^2 *
#         (λ2 - μ(α) - ν(α) + λ1*p)^2
#     )
#     return real(A - B/C)
# end

# function equilibrium_probability(λ1, λ2, α, μ, ν; tolerance = 10e-9, first_actor = 1, returnall = false)
#     if first_actor == 1
#         parameters = (λ1, λ2, α, μ, ν)
#         swapped_parameters = (λ2, λ1, α, μ, ν)
#     else
#         parameters = (λ2, λ1, α, μ, ν)
#         swapped_parameters = (λ1, λ2, α, μ, ν)
#     end
    
#     root, ismin = classify_roots(p -> diff_best_responded_cost(p, parameters...), find_valid_roots(p -> diff_best_responded_cost(p, parameters...)))
#     if sum(ismin) == 0
#         val_at_zero = analytic_best_responded_cost(0,  parameters...)
#         val_at_one = analytic_best_responded_cost(0,  parameters...)
#         root = minimum((val_at_zero, val_at_one))
#     else
#         root = first(root[ismin])
#     end
#     first_pstar(p) = pstar(p, parameters...)
#     second_pstar(p) = pstar(p, swapped_parameters...)

#     diff = 100
#     response1 = 0
#     ps = Float64[]
#     qs = Float64[]
#     if tolerance < 10e-13
#         # println("tolerance of less than 10e-13 is dangerous")
#         tolerance = 10e-14
#     end
#     while diff > tolerance
#         response1 = second_pstar(root)
#         response2 = first_pstar(response1)
#         push!(ps, response2)
#         push!(qs, response2)
#         diff = abs(response2 - root)
#         root = response2
#     end
#     if returnall == true
#         return ps, qs
#     else
#         return root, response1
#     end
# end

# Static Graph
λ1s = 0.001:0.001:1.755

probability_tuples = equilibrium_probability.(λ1s, λ2, 1, μ, ν; tolerance = 10e-9)
response_probs = first.(probability_tuples)
response_qrobs = last.(probability_tuples)

plot(λ1s, response_probs, label = "p⋆", ylims = (0,1), xlims = (0,1.76))
plot!(λ1s, response_qrobs, label = "q⋆")
plot!(xlabel = "λ1", ylabel = "probability",
    title = "optimal probabilities for various λ1")

# Animation
λ2s = 0:0.01:2
animated_λs = @animate for λ2 in λ2s
    λ1s = 0:0.01:2
    probability_tuples = equilibrium_probability.(λ1s, λ2, 1, μ, ν; tolerance = 10e-9)
    response_probs = first.(probability_tuples)
    response_qrobs = last.(probability_tuples)
    plot(title = "Traffic at each server for varying λ1 for λ2 = $(dp_string(λ2, 3))",
        legend = :outerright,
        size = (800, 600))
    plot!(λ1s, response_probs, label = "p⋆", ylims = (0,1), xlims = (0,1.76))
    plot!(λ1s, response_qrobs, label = "q⋆")
    plot!(xlabel = "λ1", ylabel = "probability")

    p1_assumptions = calculate_assumptions(response_probs, response_qrobs, λ1s, λ2, α, μ, ν; which = "both")
    p2_assumptions = calculate_assumptions(response_probs, response_qrobs, λ2, λ1s, α, μ, ν; which = "both")
    μ1_shapes = shapes_from_intervals(shaded_regions(p1_assumptions[1], λ1s))
    ν1_shapes = shapes_from_intervals(shaded_regions(p1_assumptions[2], λ1s))
    μ2_shapes = shapes_from_intervals(shaded_regions(p2_assumptions[1], λ1s))
    ν2_shapes = shapes_from_intervals(shaded_regions(p2_assumptions[2], λ1s))
    vec_shapes = [μ1_shapes, ν1_shapes, μ2_shapes, ν2_shapes]
    resized_shapes = resize_interval_shapes(vec_shapes)

    plot!(resized_shapes[1], colour=:red, label = ["μ overload p1" fill(nothing, 1, length(vec_shapes)-1)], alpha=0.3)
    plot!(resized_shapes[2], colour=:blue, label = ["ν overload p1" fill(nothing, 1, length(vec_shapes)-1)], alpha=0.3)
    plot!(resized_shapes[3], colour=:green, label = ["μ overload p2" fill(nothing, 1, length(vec_shapes)-1)], alpha=0.3)
    plot!(resized_shapes[4], colour=:yellow, label = ["ν overload p2" fill(nothing, 1, length(vec_shapes)-1)], alpha=0.3)
end

gif(animated_λs, "altering_arrival_rates.gif", fps = 30)




plot(0.5:0.0001:0.8, p -> diff_best_responded_cost(p, 1.75, 1, 1, μ, ν), ylims = (-500, 500))
plot!(0.5:0.001:0.8, p -> analytic_best_responded_cost(p, 1.75, 1, 1, μ, ν))
hline!([0], label = "", colour = :black)

begin
    roots, ismin = classify_roots(p -> diff_best_responded_cost(p, 1.75, 1, 1, μ, ν), find_valid_roots(p -> diff_best_responded_cost(big(p), 1.75, 1, 1, μ, ν)))
    vline!(roots[ismin], colour = :black, label = "min")  
    vline!(roots[.!ismin], colour = :red, label = "max")  
end

roots[ismin]

difficult_params = (1.75, 1, 1, μ, ν)
difficult_p = equilibrium_probability(difficult_params...)[1]

plot(0:0.001:1, p -> best_responded_cost(p, difficult_params...), ylims = (-50, 50))
vline!([difficult_p])

difficult_br_costs = best_responded_cost.(0.5:0.001:1, difficult_params...)
min_index = findmin(difficult_br_costs)[2]
min_at = 0.5+0.001*318

plot(0.5:0.001:1, p -> best_responded_cost(p, difficult_params...))
vline!([min_at])


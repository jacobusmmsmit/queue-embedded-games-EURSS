using ForwardDiff
using Roots
using StatsBase
using Plots
using Measures

include("helper_functions.jl")

begin
    λ1 = 1
    λ2 = 1.5
    α = 1
    β = 1
    μ(α) = 1 / (1.5 * (α))
    ν(α) = 1 / (0.7 * (α))
    parameters = (λ1, λ2, α, μ, ν)
    he_parameters = (λ1, λ2, α, μ, ν, β)
    opponent_parameters = (λ2, λ1, α, μ, ν)
    opponent_he_parameters = (λ2, λ1, α, μ, ν, β)
    qs = 0:0.001:1
end

plot(qs, he_pstar.(qs, he_parameters...))
plot!(qs, pstar.(qs, parameters...))

total_cost(he_equilibrium_probability(he_parameters...)..., parameters...)
total_cost(equilibrium_probability(parameters...)..., parameters...)

plot(qs, he_cost.(he_pstar.(qs, he_parameters...), qs, he_parameters...), ylims = (0, 50))
plot!(qs, cost.(pstar.(qs, parameters...), qs, parameters...))
plot!(qs, he_cost.(pstar.(qs, parameters...), qs, parameters..., 1))

iter = 0.001:0.001:2

begin
	probability_tuples = equilibrium_probability.(λ1, λ2, iter, μ, ν)
    he_probability_tuples = he_equilibrium_probability.(λ1, λ2, iter, μ, ν, β)
	response_probs = first.(probability_tuples)
	response_qrobs = last.(probability_tuples)
    he_response_probs = first.(he_probability_tuples)
	he_response_qrobs = last.(he_probability_tuples)

    plot(title = "Traffic at each server for varying α, λ = 1, and λ2 = $(dp_string(λ2, 3S))",
        legend = :outerright,
        size = (800, 600))
    plot!(iter, response_probs, label = "p⋆", lc = :blue, lw = 1.2)
    plot!(iter, response_qrobs, label = "q⋆", lc = :red, lw = 1.2)
    plot!(iter, he_response_probs, label = "he_p⋆", ls = :dash, lc = :blue, lw = 1.2)
    plot!(iter, he_response_qrobs, label = "he_q⋆", ls = :dash, lc = :red, lw = 1.2)
    plot!(ylims = (0,1), xlims = (0,1.76))
    plot!(xlabel = "α", ylabel = "probability")
    plot!(margin = 2mm)
end

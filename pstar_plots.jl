using Measures

include("helper_functions.jl")

begin
    λ1 = 1.4
    λ2 = 1
    α = 1
    μ(α) = 1 / (1.5 * (α))
    ν(α) = 1 / (0.7 * (α))
    parameters = (λ1, λ2, α, μ, ν)
    qs = 0:0.001:1
end

plot(qs, analytic_pstar.(qs, parameters...), label = "Best response (p*)", lw = 1.5)
plot!(xlims = (0, 1), ylims = (0, 1))
μ_boolvec, ν_boolvec = calculate_assumptions(analytic_pstar.(qs, parameters...), qs, parameters...)
plot!(shapes_from_intervals(shaded_regions(μ_boolvec, qs)), alpha = 0.3, fc = :red, label = "Assumptions violated")
plot!(xlabel = "Opponent probability (q)", ylabel = "Own probability (p)")
plot!(size=(600,300), margin = 2mm)

# savefig("best_response_plot.pdf")
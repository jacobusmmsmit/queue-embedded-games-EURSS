using Plots

include("helper_functions.jl")

begin
    λ1 = 1
    λ2 = 1
    α = 1
    μ(α) = 1 / (1.5 * (α))
    ν(α) = 1 / (0.7 * (α))
    qs = 0:0.01:1
    parameters = (λ1, λ2, α, μ, ν)
end


ps = numerical_pstar.(qs, parameters...)
a_ps = pstar.(qs, parameters...)
costs = cost.(ps, qs, parameters...)
a_costs = cost.(a_ps, qs, parameters...)

ν_assumption = ν.(α) .> λ1 .* ps .+ λ2 .* qs
μ_assumption = μ.(α) .> λ1 .* (1 .- ps)

a_ν_assumption = ν.(α) .> λ1 .* a_ps .+ λ2 .* qs
a_μ_assumption = μ.(α) .> λ1 .* (1 .- a_ps)


plot(qs, ps, label = "autodiff")
plot!(qs, a_ps, label = "analytic")

# plot!(qs, ν_assumption, colour = :red, label = "pu o ν")
# plot!(qs, μ_assumption, colour = :black, label = "pr o μ")
# plot!(qs, a_ν_assumption, colour = :red, label = "pu o ν")
# plot!(qs, a_μ_assumption, colour = :black, label = "pr o μ")

plot(qs, costs, label = "autodiff_costs")
plot!(qs, a_costs, label = "analytic_costs")
plot!(ylims = (-50, 50))

p = plot(qs, qs, (p, q) -> clamp(cost(p, q, parameters...), -50, 50), st = :surface, zlims =  (-50, 50))

# begin
#     spin_plot = @animate for i in vcat(0:0.5:90, reverse(0:0.5:90))
#         plot(p, camera = (i,i))
#     end
#     gif(spin_plot, "spin_plot.gif", fps = 60)
# end

# gif(spin_plot, "spin_plot.gif", fps = 60)

plot!(shapes_from_intervals(shaded_regions(μ_assumption, qs)), colour=:red, label = "μ violated", alpha=0.3)
plot!(shapes_from_intervals(shaded_regions(ν_assumption, qs)), colour=:blue, label = "ν violated", alpha=0.3)
plot!(shapes_from_intervals(shaded_regions(a_μ_assumption, qs)), colour=:yellow, label = "a_μ violated", alpha=0.3)
plot!(shapes_from_intervals(shaded_regions(a_ν_assumption, qs)), colour=:gray, label = "a_ν violated", alpha=0.3)


plot!(shapes_from_intervals(push!(shaded_regions(μ_assumption, qs), [0, 0.1])), colour=:red, label = "μ violated", alpha=0.3)
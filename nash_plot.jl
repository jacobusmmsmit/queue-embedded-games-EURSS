using Measures
include("helper_functions.jl")

begin
    λ1 = 1
    λ2 = 1
    α = 1
    μ(α) = 1 / (1.5 * (α))
    ν(α) = 1 / (0.7 * (α))
    parameters = (λ1, λ2, α, μ, ν)
    opponent_parameters = (λ2, λ1, α, μ, ν)
    p, q = equilibrium_probability(parameters...)
    iter = 0:0.001:1
end

cost(p, q, parameters...) ≈ cost(q, p, opponent_parameters...)
costs = cost.(iter, q, parameters...)

function find_asymptotes(vec, iter)
    sv = sign.(vec)
    iter[sv .!= circshift(sv, 1)]
end

asymptotes = find_asymptotes(costs, iter)

nash_plot = plot(iter[costs .> 0], costs[costs .> 0], label = "Response Time")
plot!(ylims = (0, 25), xlims = (0, 1), size = (600, 250), palette = :tab10)
plot!([p, p], [0, cost(p, q, parameters...)], ls = :dash, lc = :black, label = "Minimum")
plot!([0, p], [cost(p, q, parameters...), cost(p, q, parameters...)], ls = :dash, lc = :black, label = "")
vline!([asymptotes], ls = :dash, label = "Asymptotes")
plot!(xlabel = "Proportion Sent to Public Server", ylabel = "Expected Response Time", legend = :topleft, margin = 2.5mm)
savefig("outputs/nash_plot_poster.pdf")
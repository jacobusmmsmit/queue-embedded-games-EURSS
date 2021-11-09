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


itr = 0.1:0.05:2.1

poa(1.5, 1.5, α, μ, ν)
heatmap(itr, itr, (x, y) -> poa(x, y, α, μ, ν))
surface(itr, itr, (x, y) -> poa(x, y, α, μ, ν))
# heatmap(itr, itr, (x, y) -> total_cost(equilibrium_probability(x, y, α, μ, ν)..., x, y, α, μ, ν))
# heatmap(itr, itr, (x, y) -> total_cost(he_equilibrium_probability(x, y, α, μ, ν, 1)..., x, y, α, μ, ν))
# plot(0.666:0.00001:0.667, 2.095:0.00001:2.0953, (x, y) -> poa(x, y, α, μ, ν), st = :contourf)
poa(0.66655, 2.09522, 1, μ, ν)
poa(0.666, 2.095, 1, μ, ν)
poa(0.66, 2.0, 1, μ, ν)
poac(x, y) = poa(x, y, 1, μ, ν)
poac(0.1, 2.095)

plot!(
    xlabel = "Player 1 arrival rate",
    ylabel = "Player 2 arrival rate",
    title = "Price of Anarchy for Various Arrival Rates",
    xticks = 0:0.25:2,
    yticks = 0:0.25:2,
    clim = (1, Inf),
    margin = 6pt,
    size = (450, 400),
)
annotate!(1.45, 1.45, text("Infeasible \n Arrival Rates", 9, rotation = -45))

savefig("outputs/price_of_anarchy_heatmap.pdf")
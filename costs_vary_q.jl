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
    qs = 0:0.001:1
end

p1_cost = cost.(pstar.(qs, parameters...), qs, parameters...)
p2_cost = cost.(qs, pstar.(qs, parameters...), opponent_parameters...)
p2_mask = p2_cost .>= 0

plot(qs, p1_cost, label = "Player 1 cost")
plot!(qs[p2_mask], p2_cost[p2_mask], label = "Player 2 cost")
vline!([p], label = "Equilibrium probability", lc = :red, ls = :dash)
xlabel!("Opponent's probability")
ylabel!("Cost")
plot!(size = (600, 250), ylims = (0, 25), xlims = (0.2,1), margin = 2.5mm)

savefig("outputs/costs_vary_q.pdf")
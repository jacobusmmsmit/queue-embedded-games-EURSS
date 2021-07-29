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

p1_cost = cost.(analytic_pstar.(qs, parameters...), qs, parameters...)
p2_cost = cost.(qs, analytic_pstar.(qs, parameters...), opponent_parameters...)
p2_mask = p2_cost .>= 0
p1_mask = p1_cost .>= 0

begin
    cost_plot = plot(qs[p1_mask], p1_cost[p1_mask], label = "Player 1 cost")
    plot!(qs[p2_mask], p2_cost[p2_mask], label = "Player 2 cost")
    vline!([p], label = "Equilibrium probability", lc = :red, ls = :dash)
    xlabel!("Opponent's probability")
    ylabel!("Cost")
    plot!(size = (600, 250), ylims = (0, 25), xlims = (0.2,1), margin = 2.5mm)
end



# savefig("outputs/costs_vary_q.pdf")
q
perturbed_q = 0.6
current_cost = cost(p, q, parameters...)
p1_newcost = cost(analytic_pstar(perturbed_q, parameters...), perturbed_q, parameters...)
analytic_best_responded_cost(q, parameters...)
p2_newcost = cost(perturbed_q, analytic_pstar(perturbed_q, parameters...), opponent_parameters...)

r2perturb = analytic_pstar(perturbed_q, parameters...)
r3p = analytic_pstar(r2perturb, opponent_parameters...)
r4p = analytic_pstar(r3p, parameters...)

ps = analytic_pstar.(qs, parameters...)

p1_plot = plot(qs, ps, label = "p*", lims = (0, 1))
vline!(p1_plot, [perturbed_q], ls = :dash)
vline!(p1_plot, [r2perturb], ls = :dash)
vline!(p1_plot, [r3p], ls = :dash)
vline!(p1_plot, [r4p], ls = :dash)

p1_plot = plot(qs, analytic_pstar.(qs, opponent_parameters...), label = "q*", lims = (0, 1))
vline!([q], label = "Pareto probability", lc = :red, ls = :dash)
xlabel!("Opponent's probability")
ylabel!("Cost")

p1_queue_length = (1 .- ps) .* λ1 .* private_cost.(ps, qs, parameters...)
p2_queue_length = (1 .- qs) .* λ2 .* private_cost.(qs, ps, opponent_parameters...)
public_queue_length = ((λ1 .* ps) + (λ2 .* qs)) .* public_cost.(ps, qs, parameters...)

plot(qs, p1_queue_length, label = "Player 1 queue length")
plot!(qs[p2_queue_length .> 0], p2_queue_length[p2_queue_length .> 0], label = "Player 2 queue length")
plot!(qs, public_queue_length, label = "Public queue length")
ylims!(0,30)


q
analytic_pstar(q, parameters...)


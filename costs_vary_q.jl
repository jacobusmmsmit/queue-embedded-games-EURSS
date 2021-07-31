using Measures
include("helper_functions.jl")

begin
    λ1 = 1
    λ2 = 1.5
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
p1_mask = p1_cost .>= 0

begin
    cost_plot = plot(qs[p1_mask], p1_cost[p1_mask], label = "Player 1 cost")
    plot!(qs[p2_mask], p2_cost[p2_mask], label = "Player 2 cost")
    vline!([q], label = "Equilibrium", lc = :red, ls = :dash)
    xlabel!("Player 2 probability")
    ylabel!("Cost")
    plot!(size = (600, 250), ylims = (5, 20), xlims = (0.57, 0.68), margin = 2.5mm)
    plot!(legend = :bottomright)
end
savefig("outputs/costs_vary_q_15.pdf")

begin
    total_cost_var = p1_cost + p2_cost
    p1gp2 = sign.(p1_cost[p1_mask] - p2_cost[p1_mask])
    tc_plot = plot(cost_plot, qs, total_cost_var, label = "Total cost")
    plot!(xlims = (0.57, 0.65), ylims = (0, 40))
    mincost, pos = findmin(total_cost_var[0.56 .< qs .< 0.7])
    vline!([qs[0.56 .< qs .< 0.7][pos]], label = "Lower total cost", ls = :dash, lc = :black)
    vline!([qs[p1_mask][p1gp2 .!= circshift(p1gp2, 1)]], label = "Minmax cost", ls = :dash)
    plot!(size = (600, 400), legend = :topright)
end
savefig("outputs/soc_opt_15.pdf")



calculate_assumptions(0.5, 0.5, parameters...)

# function total_cost(p, q, λ1, λ2, α, μ, ν)
#     p1_cost = cost(p, q, λ1, λ2, α, μ, ν)
#     p2_cost = cost(q, p, λ2, λ1, α, μ, ν)
    
#     return (p1_cost, p2_cost, (p1_cost > 0 ? p1_cost : 100000) + (p2_cost > 0 ? p2_cost : 100000))
# end



# begin
#     qn = q - 0.004
#     pn = pstar(qn, parameters...)
#     total_cost(p, q, parameters...) .- total_cost(pn, qn, parameters...)
#     # println("q: ", q)
#     # println("qn: ", qn)
#     # println("p: ", p)
#     # println("pn: ", pn)
# end

total_cost(p, q, parameters...)
total_cost(p, 0.5, parameters...)
total_cost(pn, qn, parameters...)

res = total_cost.(qs, pstar.(qs, parameters...), parameters...)
p1_c = getindex.(res, 1)
p2_c = getindex.(res, 2)
tcs = getindex.(res, 3)
plot(qs, p1_c)
plot!(qs, p2_c)
plot!(qs, tcs)
plot!(ylims = (-50, 50))
vline!([pn])
vline!([q])


plot(
    0:0.005:1,
0:0.005:1,
(p, q) -> total_cost(p, q, parameters...),
    st = :surface,
    zlims = (-50, 50),
    camera = (90, 90))
scatter!([p], [q], [total_cost(p, q, parameters...)], legend = nothing)

2
# q
# perturbed_q = 0.6
# current_cost = cost(p, q, parameters...)
# p1_newcost = cost(pstar(perturbed_q, parameters...), perturbed_q, parameters...)
# analytic_best_responded_cost(q, parameters...)
# p2_newcost = cost(perturbed_q, pstar(perturbed_q, parameters...), opponent_parameters...)

# r2perturb = pstar(perturbed_q, parameters...)
# r3p = pstar(r2perturb, opponent_parameters...)
# r4p = pstar(r3p, parameters...)

# ps = pstar.(qs, parameters...)

# p1_plot = plot(qs, ps, label = "p*", lims = (0, 1))
# vline!(p1_plot, [perturbed_q], ls = :dash)
# vline!(p1_plot, [r2perturb], ls = :dash)
# vline!(p1_plot, [r3p], ls = :dash)
# vline!(p1_plot, [r4p], ls = :dash)

# p1_plot = plot(qs, pstar.(qs, opponent_parameters...), label = "q*", lims = (0, 1))
# vline!([q], label = "Pareto probability", lc = :red, ls = :dash)
# xlabel!("Opponent's probability")
# ylabel!("Cost")

# p1_queue_length = (1 .- ps) .* λ1 .* private_cost.(ps, qs, parameters...)
# p2_queue_length = (1 .- qs) .* λ2 .* private_cost.(qs, ps, opponent_parameters...)
# public_queue_length = ((λ1 .* ps) + (λ2 .* qs)) .* public_cost.(ps, qs, parameters...)

# plot(qs, p1_queue_length, label = "Player 1 queue length")
# plot!(qs[p2_queue_length .> 0], p2_queue_length[p2_queue_length .> 0], label = "Player 2 queue length")
# plot!(qs, public_queue_length, label = "Public queue length")
# ylims!(0,30)


# q
# pstar(q, parameters...)


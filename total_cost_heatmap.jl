using ForwardDiff
using ReverseDiff
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
    step = 0.0005
    xs = 0.3:step:0.65
    ys = 0.5:step:0.75
end


threeD_cost = (p, q) -> total_cost(p, q, parameters...)
heatcost = plot(xs, ys, threeD_cost, st = :heatmap, clims = (0, 150))
plot!(xlims = (0.3, 0.6), ylims = (0.5, 0.75))
plot!([p], [q], st = :scatter, shape = :circle, msc = :black, label = "Equilibrium")
plot!(xlabel = "Player 1 probability", ylabel = "Player 2 probability")

begin
    gtape = ReverseDiff.GradientTape(vec -> threeD_cost(vec...), [p, q]) 
    out = zeros(2)
    dtotal(ϕ, ψ) = ReverseDiff.gradient!(out, gtape, [ϕ, ψ])
end

difference =  1
p1, q1 = p, q
while difference > 1e-6
    p0, q0 = p1, q1
    proots = first(find_valid_roots(a -> first(dtotal(a, q0))))
    p1 = proots[findmin(abs.(proots .- p0))[2]]
    qroots = first(find_valid_roots(b -> last(dtotal(p1, b))))
    q1 = qroots[findmin(abs.(qroots .- q0))[2]] 
    difference = maximum([abs(p1 - p0), abs(q1 - q0)])
end

plot!(heatcost, [p1], [q1], st = :scatter, label = "Minimum Total Cost")

begin
    possibleroots = first(find_valid_roots(q -> ForwardDiff.derivative(ψ -> total_cost(pstar(ψ, parameters...), ψ, parameters...), q)))
    roots = possibleroots[ForwardDiff.derivative.(q -> ForwardDiff.derivative(ψ -> total_cost(pstar(ψ, parameters...), ψ, parameters...), q), possibleroots) .> 0]
    q_p1br = roots[findmin(abs.(roots .- q))[2]]
    p_p1br = pstar(q_p1br, parameters...)
    plot!(heatcost, [p_p1br], [q_p1br], st = :scatter, label = "Min where P1 best responds")
end

begin
    possibleroots = first(find_valid_roots(p -> ForwardDiff.derivative(ψ -> total_cost(ψ, qstar(ψ, parameters...), parameters...), p)))
    roots = possibleroots[ForwardDiff.derivative.(p -> ForwardDiff.derivative(ψ -> total_cost(ψ, qstar(ψ, parameters...), parameters...), p), possibleroots) .> 0]
    p_p2br = roots[findmin(abs.(roots .- p))[2]]
    q_p2br = qstar(p_p2br, parameters...)
    plot!(heatcost, [p_p2br], [q_p2br], st = :scatter, label = "Min where P2 best responds")
end
plot!(colorbar_title = "Total Cost")
savefig("outputs/total_cost_heatmap.pdf")

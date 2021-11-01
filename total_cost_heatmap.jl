using ForwardDiff
using ReverseDiff
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
    step = 0.0005
    xs = 0.3:step:0.65
    ys = 0.5:step:0.75
end


threeD_cost = (p, q) -> total_cost(p, q, parameters...)
heatcost = plot(xs, ys, threeD_cost, st = :contourf)
plot!(xlims = (0.4, 0.45), ylims = (0.58, 0.62), clims = (20, 27))
scattersize = 5
plot!([p], [q], st = :scatter, shape = :circle, msc = :black, ms = scattersize, label = "Equilibrium")
plot!(xlabel = "Agent 1 proportion", ylabel = "Agent 2 proportion")

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

plot!(heatcost, [p1], [q1], st = :scatter, ms = scattersize, label = "Minimum Response Time")
plot!(colorbar_title = "Sum of Average Response Time")
plot!(size = (600, 300), margin = 1.5mm)

heatcost
begin
    possibleroots = first(find_valid_roots(q -> ForwardDiff.derivative(ψ -> total_cost(pstar(ψ, parameters...), ψ, parameters...), q)))
    roots = possibleroots[ForwardDiff.derivative.(q -> ForwardDiff.derivative(ψ -> total_cost(pstar(ψ, parameters...), ψ, parameters...), q), possibleroots) .> 0]
    q_p1br = roots[findmin(abs.(roots .- q))[2]]
    p_p1br = pstar(q_p1br, parameters...)
    heatcost = plot(heatcost, [p_p1br], [q_p1br], st = :scatter, ms = scattersize, label = "Min where A1 best responds")
end

begin
    possibleroots = first(find_valid_roots(p -> ForwardDiff.derivative(ψ -> total_cost(ψ, qstar(ψ, parameters...), parameters...), p)))
    roots = possibleroots[ForwardDiff.derivative.(p -> ForwardDiff.derivative(ψ -> total_cost(ψ, qstar(ψ, parameters...), parameters...), p), possibleroots) .> 0]
    p_p2br = roots[findmin(abs.(roots .- p))[2]]
    q_p2br = qstar(p_p2br, parameters...)
    plot!(heatcost, [p_p2br], [q_p2br], st = :scatter, ms = scattersize, label = "Min where A2 best responds")
end

savefig("outputs/total_cost_contour_poster.pdf")

probability_path = he_equilibrium_probability.(λ1, λ2, α, μ, ν, 0:0.001:1)

begin
    heatcost_path = plot(heatcost, legend = :topright, size = (500, 300), xlims = (0.41, 0.44), ylims = (0.59, 0.615))
    n = 5
    for i in 1:n
        interval = Int(floor(length(probability_path) / n))
        slice = 1 + (i - 1) * interval:i * interval
        plot!(heatcost_path,
            first.(probability_path[slice]),
            last.(probability_path[slice]),
            lc = :black,
            lw = 1,
            ls = :solid,
            arrow = i == n ? false : true,
            label = i == 1 ? "Homo-Egualis path" : "")
        # if i != n
        #     annotate!(
        #         heatcost_path,
        #         probability_path[i*interval + 1]...,
        #         text(
        #             round(collect(0:0.001:1)[i*interval], digits = 2),
        #             :bottom,
        #             :green,
        #             11)
        #         )
        # end
    end
    heatcost_path
end

savefig("outputs/homo-egualis_path_contour_poster.pdf")
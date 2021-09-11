using ForwardDiff
using ReverseDiff
include("helper_functions.jl")

begin
    λ1 = 0.4
    λ2 = 1.5
    α = 1
    μ(α) = 1 / (1.5 * (α))
    ν(α) = 1 / (0.7 * (α))
    parameters = (λ1, λ2, α, μ, ν)
    opponent_parameters = (λ2, λ1, α, μ, ν)
    p, q = equilibrium_probability(parameters...)
    step = 0.0005
    xs = 0:step:1
    ys = 0:step:1
end

begin
    threeD_cost = (p, q) -> total_cost(p, q, parameters...)
    heatcost = plot(xs, ys, threeD_cost, st = :contourf)
    scattersize = 5
    plot!([p], [q], st = :scatter, shape = :circle, msc = :black, ms = scattersize, label = "Equilibrium")
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

    multiplier = 0.1
    plot_xlims = p1 > p ? ((1-multiplier) * p, (1+multiplier) * p1) : ((1-multiplier) * p1, (1+multiplier) * p) #(0.4, 0.45)
    plot_ylims = q1 > q ? ((1-multiplier) * q, (1+multiplier) * q1) : ((1-multiplier) * q1, (1+multiplier) * q) #(0.58, 0.62)
    plot_clims = (threeD_cost(p1, q1)*(1-multiplier), threeD_cost(p, q)*(1+multiplier))

    plot!(heatcost, [p1], [q1], st = :scatter, ms = scattersize, label = "Minimum Total Cost")
    plot!(xlims = plot_xlims, ylims = plot_ylims, clims = plot_clims)
    plot!(colorbar_title = "Total Cost")
    plot!(size = (500, 300))
end

begin
    possibleroots = first(find_valid_roots(q -> ForwardDiff.derivative(ψ -> total_cost(pstar(ψ, parameters...), ψ, parameters...), q)))
    roots = possibleroots[ForwardDiff.derivative.(q -> ForwardDiff.derivative(ψ -> total_cost(pstar(ψ, parameters...), ψ, parameters...), q), possibleroots) .> 0]
    q_p1br = roots[findmin(abs.(roots .- q))[2]]
    p_p1br = pstar(q_p1br, parameters...)
    heatcost_extrapoints = plot(heatcost, [p_p1br], [q_p1br], st = :scatter, ms = scattersize, label = "Min where P1 best responds")
end

begin
    possibleroots = first(find_valid_roots(p -> ForwardDiff.derivative(ψ -> total_cost(ψ, qstar(ψ, parameters...), parameters...), p)))
    roots = possibleroots[ForwardDiff.derivative.(p -> ForwardDiff.derivative(ψ -> total_cost(ψ, qstar(ψ, parameters...), parameters...), p), possibleroots) .> 0]
    p_p2br = roots[findmin(abs.(roots .- p))[2]]
    q_p2br = qstar(p_p2br, parameters...)
    plot!(heatcost_extrapoints, [p_p2br], [q_p2br], st = :scatter, ms = scattersize, label = "Min where P2 best responds")
end



# savefig("outputs/total_cost_contour.pdf")

probability_path = he_equilibrium_probability.(λ1, λ2, α, μ, ν, 0:0.001:1)

begin
    heatcost_path = plot(heatcost, legend = :topright, size = (500, 300), xlims = plot_xlims, ylims = plot_ylims)
    n = 5
    for i in 1:n
        interval = Int(floor(length(probability_path)/n))
        slice = 1 + (i-1) * interval: i*interval
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

# savefig("outputs/homo-egualis_path_contour.pdf")
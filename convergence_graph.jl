using Measures
include("helper_functions.jl")

begin
    λ1 = 0.4
    λ2 = 1.4
    α = 1
    μ(α) = 1 / (1.5 * (α))
    ν(α) = 1 / (0.7 * (α))
    parameters = (λ1, λ2, α, μ, ν)
end

ps, qs = equilibrium_probability(parameters...; returnall=true)

cost(plim, qlim, parameters...)

plim, qlim = equilibrium_probability(parameters...; returnall=false)

ps_plot = plot(ylabel = "Probability", xlabel = "", size = (600,250), legend = :topright,lw = 1.5)
plot!(ps_plot, ps, label = "Best response", lc = :blue)
hline!([plim], lc = :red, ls = :dash, label = "Optimal probability")

qs_plot = plot(qs, label = "Opponent best response", size = (600,250), legend = :bottomright, lc = :orange, lw = 1.5, ylabel = "Probability")
hline!([qlim], lc = :red, ls = :dash, label = "Optimal probability")

plot(ps_plot, qs_plot, xlabel = "Number of Iterations", margin = 2mm)
savefig("convergence_plot.pdf")
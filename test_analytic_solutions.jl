using Plots
using ForwardDiff
using Roots
include("helper_functions.jl")

begin
    λ1 = 0.5
    λ2 = 1
    α = 1
    μ(α) = 1 / (1.5 * (α))
    ν(α) = 1 / (0.7 * (α))
    parameters = (λ1, λ2, α, μ, ν)
    q = 0.5
    qs = 0:0.001:1
end

# Test 1: Do the derivatives agree?
# function adcost(p, q, λ1, λ2, α, μ, ν)
#     (λ1 * p) / (λ1 * p - ν(α) + λ2 * q)^2 -
#     1/(μ(α) - λ1 + λ1 * p) - 
#     1/(λ1 * p - ν(α) + λ2 * q) +
#     (λ1 * (p - 1))/(μ(α) - λ1 + λ1 * p)^2
# end

# dcost = ∂cost_∂p

ps = 0:0.001:1
plot(ps, ∂cost_∂p.(ps, q, parameters...))
plot!(ps, diff_cost.(ps, q, parameters...))
plot!(ylims = (-500, 500))

# @benchmark adcost(0.5, q, parameters...)
# @benchmark dcost(0.5, q, parameters...)

# Test 2: Do the points where derivative is zero agree? (aka pstar)
sol, _ = first.(find_valid_roots(p -> ∂cost_∂p(p, q, parameters...)))
# function apstar(q, λ1, λ2, α, μ, ν)
#     p = (λ1*Complex(μ(α)*(ν(α) - λ2*q))^(1/2) - μ(α)*Complex(μ(α)*(ν(α) - λ2*q))^(1/2) - ν(α)*Complex(μ(α)*(ν(α) - λ2*q))^(1/2) - λ1*ν(α) + 2*μ(α)*ν(α) + λ2*q*Complex(μ(α)*(ν(α) - λ2*q))^(1/2) + λ1*λ2*q - 2*λ2*μ(α)*q)/(λ1*μ(α) - λ1*ν(α) + λ1*λ2*q) |> real
#     return clamp(p, 0, 1)
# end
asol = analytic_pstar(q, parameters...)

diff_cost(sol, q, parameters...)
diff_cost(asol, q, parameters...)

# Test 3: Do the pstars agree?
# function pstar(q, λ1, λ2, α, μ, ν)
#     roots = find_valid_roots(p -> ∂cost_∂p(p, q, λ1, λ2, α, μ, ν))
#     if length(roots[1]) > 1
#         root_values, root_is_minimum =
#             classify_roots(p -> ∂cost_∂p(p, q, λ1, λ2, α, μ, ν), roots)
#         return (root_values[root_is_minimum])[1]
#     else
#         return (roots[1])[1]
#     end
# end


qs = 0:0.0001:1
ps = pstar.(qs, parameters...)
a_ps = analytic_pstar.(qs, parameters...)

response_plot = plot(qs, ps, label = "autodiff")
plot!(qs, a_ps, label = "analytic")
plot_assumptions!(response_plot, [ps, a_ps], qs, parameters...)
plot!(response_plot, legend = :outertopleft, size = (690, 400))

costs = cost.(ps, qs, parameters...)
a_costs = cost.(a_ps, qs, parameters...)

cost_plot = plot(qs, costs, label = "autodiff")
plot!(qs, a_costs, label = "analytic")
plot!(ylims = (-50, 50))
plot_assumptions!(cost_plot, [ps, a_ps], qs, parameters..., height = 10)

traffic_plot = plot(
    title = "Traffic at each server for varying q",
    xlabel = "Opponent probability (q)",
    ylabel = "Intensity of arrivals",
    legend = :topleft,
    xlims = (0, 1),
    ylims = ylims = (0, Inf)
    )

plot!(qs, λ1 .* a_ps .+ λ2 .* qs, label = "ν(α)")
hline!([ν.(α)], label = "ν(α) limit", ls = :dash, lw = 1.5)

plot!(qs, λ1 .* (1 .- a_ps), label = "μ(α)")
hline!([μ(α)], label = "μ(α) limit", ls = :dash, lw = 1.5)

intersection_index = last(findmax(λ1 .* a_ps .+ λ2 .* qs .> ν(α)))
x_intersection = qs[intersection_index]
y_intersection_max = maximum([
    (λ1 .* (1 .- a_ps))[intersection_index],
    (λ1 .* a_ps .+ λ2 .* qs)[intersection_index]]
    )
plot!([x_intersection, x_intersection], [0, y_intersection_max], lw = 1, lc = :maroon, label = "")

traffic_plot_anim = @animate for α ∈ vcat(collect(0.5:0.005:2.1), collect(reverse(0.5:0.005:2.1)))
    a_ps = apstar.(qs, (λ1, λ2, α, μ, ν)...)
    traffic_plot = plot(
    title = "Traffic at each server for varying qs at α = $(dp_string(α, 3))",
    # xlabel = "Opponent probability (q)",
    ylabel = "Intensity of arrivals",
    legend = :topleft,
    xlims = (0, 1),
    ylims = ylims = (0, Inf),
    xticks = 0:0.2:1,
    )

    plot!(qs, λ1 .* a_ps .+ λ2 .* qs, label = "ν(α)")
    hline!([ν.(α)], label = "ν(α) limit", ls = :dash, lw = 1.5)

    plot!(qs, λ1 .* (1 .- a_ps), label = "μ(α)")
    hline!([μ(α)], label = "μ(α) limit", ls = :dash, lw = 1.5)
    
    does_intersect, intersection_index = findmax(λ1 .* a_ps .+ λ2 .* qs .> ν(α))
    x_intersection = qs[intersection_index]
    y_intersection_max = maximum([
        (λ1 .* (1 .- a_ps))[intersection_index],
        (λ1 .* a_ps .+ λ2 .* qs)[intersection_index]]
        )
    plot!([x_intersection, x_intersection], [0, y_intersection_max], ls = :dash, lw = 1, lc = :black, label = "")
    

    pstar_plot = plot(
        qs,
        a_ps,
        lc = :black,
        label = "Best Response (p⋆)",
        xlabel = "Opponent probability (q)",
        xticks = 0:0.2:1,
        ylabel = "Response probability (p)",
        xlims = (0, 1),
        ylims = (0, 1))
    plot!(
        [x_intersection, x_intersection],
        [0, a_ps[intersection_index]],
        label = "",
        ls = :dash,
        lc = :black)

    if does_intersect
        plot!(
            qs[1:intersection_index],
            zeros(length(qs[1:intersection_index])),
            fillrange = a_ps[1:intersection_index],
            fc = :green,
            fa = 0.35,
            label = "BR Feasibile")
        plot!(
            qs[intersection_index + 1:end],
            zeros(length(qs[intersection_index + 1:end])),
            fillrange = a_ps[intersection_index + 1:end],
            fc = :red,
            fa = 0.35,
            label = "BR Infeasibile")
    else
        plot!(
            qs,
            zeros(length(qs)),
            fillrange = a_ps,
            fc = :green,
            fa = 0.35,
            label = "BR Feasibile")
        # Have to put something here so it doesn't change the legend
        plot!(
            [1, 2],
            zeros(2),
            fillrange = [2, 3],
            fc = :red,
            fa = 0.35,
            label = "BR Infeasibile")
    end

    plot(traffic_plot, pstar_plot, layout = (2, 1), size = (800, 800))
end

gif(traffic_plot_anim, "traffic_plot_anim.gif", fps = 30)
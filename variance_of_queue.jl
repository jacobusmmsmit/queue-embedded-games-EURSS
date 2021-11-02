using Measures
include("helper_functions.jl")

begin
    αs = 0.2:0.01:1.57
    total_volume = 1
    ARs = total_volume ./ αs
    λ1 = 1
    λ2 = 1
    α = 1
    μ(α) = 1 / (1.5 * (α))
    ν(α) = 1 / (0.7 * (α))
    parameters = (λ1, λ2, α, μ, ν)
end

function variance_of_cost(p, q, λ1, λ2, α, μ, ν)
    (p * public_cost(p, q, λ1, λ2, α, μ, ν))^2 + ((1 - p) * private_cost(p, λ1, α, μ))^2
end

plot(αs,
variance_of_cost.(0.57, 0.57, ARs, ARs, αs, μ, ν),
    yaxis = :log10,
    yticks = 10.0.^[-1, 0, 1, 2, 3],
    legend = false,
    title = "Variance by Job Size",
    xlabel = "Job Size",
    ylabel = "Variance of Cost",
    size = [400, 300],)

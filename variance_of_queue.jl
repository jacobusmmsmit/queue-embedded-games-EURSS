using Measures
include("helper_functions.jl")

begin
    λ1 = 0.4
    λ2 = 1.4
    α = 1
    μ(α) = 1 / (1.5 * (α))
    ν(α) = 1 / (0.7 * (α))
    parameters = (λ1, λ2, α, μ, ν)
    αs = 0.5:0.01:1.5
end

function variance_of_cost(p, q, λ1, λ2, α, μ, ν)
    (p*public_cost(p, q, λ1, λ2, α, μ, ν))^2 + ((1-p) * private_cost(p, λ1, α, μ))^2
end

plot(αs, variance_of_cost.(0.5, 0.5, λ1, λ2, αs, μ, ν))

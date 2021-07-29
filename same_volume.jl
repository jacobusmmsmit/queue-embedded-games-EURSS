include("helper_functions.jl")

begin
    λ1 = 1
    λ2 = 
    α = 1
    μ(α) = 1 / (1.5 * (α))
    ν(α) = 1 / (0.7 * (α))
    parameters = (λ1, λ2, α, μ, ν)
    opponent_parameters = (λ2, λ1, α, μ, ν)
    p, q = equilibrium_probability(parameters...)
    qs = 0:0.001:1
end
αs = collect(0.5:0.01:3)
λ1 = 0.7 ./ αs
λ2 = 1.2 ./ αs

plot(αs, first.(equilibrium_probability.(λ1, λ2, αs, μ, ν)))
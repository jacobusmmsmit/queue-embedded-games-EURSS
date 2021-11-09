using Random
using StatsBase
using Distributions
using StatsPlots

include("helper_functions.jl")

begin
    λ1 = 1
    λ2 = 1.5
    α = 1
    β = 1
    μ(α) = 1 / (1.5 * (α))
    ν(α) = 1 / (0.7 * (α))
    parameters = (λ1, λ2, α, μ, ν)
    he_parameters = (λ1, λ2, α, μ, ν, β)
    # opponent_parameters = (λ2, λ1, α, μ, ν)
    # opponent_he_parameters = (λ2, λ1, α, μ, ν, β)
    # qs = 0:0.001:1
end

Hyperexponential(parameters, weights) = MixtureModel(Exponential, parameters, weights)
RTDist(p, q, λ1, λ2, α, μ, ν) = Hyperexponential([(1 / (ν(α) - λ1 * p - λ2 * q)), (1 / (μ(α) - λ1 * (1 - p)))], [p, 1 - p])
RT = RTDist(0.57, 0.57, 1, 1, 1, μ, ν)
RT2 = RTDist(0.57, 0.57, 0.5, 0.5, 2, μ, ν)
density(log.(rand(RT, 1000000)), lw = 4, label = "Theoretical Response Distribution")
density!(log.(rand(RT2, 1000000)), lw = 4, label = "Theoretical Response Distribution")

JSs = 0.1:0.01:1
AR_JS_ratio = 1
ARs = AR_JS_ratio ./ JSs 

begin
RTs = RTDist.(0.57, 0.57, ARs, ARs, JSs, μ, ν)
plot_mean = plot(JSs,
    mean.(RTs),
    title = "Mean Response Time by Job Size",
    xlabel = "Job Size",
    ylabel = "Mean Response Time",
    # label = "Theoretical Mean",
    size = [400, 300],
    legend = :none)

plot_variance = plot(JSs,
    var.(RTs),
    title = "Variance by Job Size",
    xlabel = "Job Size",
    ylabel = "Variance of Response Time",
    # label = "Theoretical Variance",
    size = [400, 300],
    legend = :none)

mean_var_plot = plot(plot_mean, plot_variance, size = [800, 300], margin = 12pt)
end

savefig("outputs/theoretical_mean_var.pdf")
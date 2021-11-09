using ForwardDiff
using Roots
using StatsBase
using Plots

function private_cost(p, λ1, α, μ)
    1 / (μ(α) - λ1 + (λ1 * p))
end

private_cost(p, q, λ1, λ2, α, μ, ν) = private_cost(p, λ1, α, μ)

"""
"""
function public_cost(p, q, λ1, λ2, α, μ, ν)
    (1 / (ν(α) - λ1 * p - λ2 * q))
end


"""
    cost(p, q, λ1, λ2, α, μ, ν)
Input parameters and return the cost
"""
function cost(p, q, λ1, λ2, α, μ, ν)
    p * public_cost(p, q, λ1, λ2, α, μ, ν) + (1-p) * private_cost(p, λ1, α, μ)
end

total_cost(p, q, λ1, λ2, α, μ, ν) = cost(p, q, λ1, λ2, α, μ, ν) + cost(q, p, λ2, λ1, α, μ, ν)

"""
    ∂cost_∂p(p, q, λ1, λ2, α, μ, ν)
Input parameters and return the derivative of the cost with respect to p at the
specified inputs. Derivative calculated using ForwardDiff, provides the same
output as `diff_cost` but is less performant.
"""
function ∂cost_∂p(p, q, λ1, λ2, α, μ, ν)
    ForwardDiff.derivative(ϕ -> cost(ϕ, q, λ1, λ2, α, μ, ν), p)
end

"""
    diff_cost(p, q, λ1, λ2, α, μ, ν)
Input parameters and return the derivative of the cost with respect to p at the
specified inputs. Derivative calculated using ForwardDiff, provides the same
output as `∂cost_∂p` but is more performant.
"""
function diff_cost(p, q, λ1, λ2, α, μ, ν)
    diff_public = (ν(α) - λ2*q)/(λ1*p - ν(α) + λ2*q)^2
    diff_private = μ(α)/(μ(α) - λ1 + λ1*p)^2
    return diff_public - diff_private
end


"""
    find_valid_roots(func)
Input a function and this will find values between 0 and 1 that minimise or
maximise that function. Both interior and boundary solutions are found. If
a boundary solution is found the second part of the return argument is true
and false if not.
"""
function find_valid_roots(func)
    roots = find_zeros(func, 0, 1)
    if length(roots) == 0
        return argmin([func(0), func(1)]) .- 1, true
    else
        return roots, false
    end
end


"""
    classify_roots(func)
Input a function (which should be a derivative) and some values at which its
output is 0 and this function returns the roots and a boolean array in which
1 means minimum and 0 means maximum.
"""
function classify_roots(func, roots)
    # roots is larger than 1 if it is a direct input from `find_valid_roots`
    # in which case it contains whether the root is interior or exterior
    if length(roots) > 1  
        if roots[2] == true
            return [roots[1]], [true]
        else
            root_values = roots[1]
        end
    else
        root_values = roots
    end
    second_derivatives = ForwardDiff.derivative.(func, root_values)
    classifications = second_derivatives .> 0 # 1 is minimum
    return root_values, classifications
end


"""
    function numerical_pstar(q, λ1, λ2, α, μ, ν)
Returns the best response to an opponents move from the perspective of the
player responding i.e. the argument λ1  may be λ2 in practice if player 2 is
responding to player 1. Calculated numerically, less accurate than `pstar`
"""
function numerical_pstar(q, λ1, λ2, α, μ, ν)
    roots = find_valid_roots(p -> diff_cost(p, q, λ1, λ2, α, μ, ν))
    if length(roots[1]) > 1
        root_values, root_is_minimum =
            classify_roots(p -> diff_cost(p, q, λ1, λ2, α, μ, ν), roots)
        return (root_values[root_is_minimum])[1]
    else
        return (roots[1])[1]
    end
end


"""
    pstar(q, λ1, λ2, α, μ, ν)
Returns the best response to an opponents move from the perspective of the
player responding i.e. the argument λ1  may be λ2 in practice if player 2 is
responding to player 1. Calculated analytically, more accurate than `numerical_pstar`
"""
function pstar(q, λ1, λ2, α, μ, ν)
    p = (λ1*Complex(μ(α)*(ν(α) - λ2*q))^(1/2) - μ(α)*Complex(μ(α)*(ν(α) - λ2*q))^(1/2) - ν(α)*Complex(μ(α)*(ν(α) - λ2*q))^(1/2) - λ1*ν(α) + 2*μ(α)*ν(α) + λ2*q*Complex(μ(α)*(ν(α) - λ2*q))^(1/2) + λ1*λ2*q - 2*λ2*μ(α)*q)/(λ1*μ(α) - λ1*ν(α) + λ1*λ2*q) |> real
    return clamp(p, 0, 1)
end

qstar(p, λ1, λ2, α, μ, ν) = pstar(p, λ2, λ1, α, μ, ν)

# function pstar(q, λ1, λ2, α, μ, ν)
#     p = (λ1*Complex(μ(α)*(ν(α) - λ2*q))^(1/2) - μ(α)*Complex(μ(α)*(ν(α) - λ2*q))^(1/2) - ν(α)*Complex(μ(α)*(ν(α) - λ2*q))^(1/2) - λ1*ν(α) + 2*μ(α)*ν(α) + λ2*q*Complex(μ(α)*(ν(α) - λ2*q))^(1/2) + λ1*λ2*q - 2*λ2*μ(α)*q)/(λ1*μ(α) - λ1*ν(α) + λ1*λ2*q)
#     if isreal(p) && (0 < Real(p) < 1)
#         return Real(p)
#     else
#         val_at_zero = analytic_best_responded_cost(0,  parameters...)
#         val_at_one = analytic_best_responded_cost(1,  parameters...)
#         return findmin((val_at_zero, val_at_one))[2]
#     end
# end


"""
    best_responded_cost(p, λ1, λ2, α, μ, ν)
Returns the cost given an input p assuming that the opponent will best respond
to this p. Calculated analytically, but less performant than
`analytic_best_responded_cost`.
"""
function best_responded_cost(p, λ1, λ2, α, μ, ν)
    cost(p, pstar(p, λ2, λ1, α, μ, ν), λ1, λ2, α, μ, ν)
end

"""
    analytic_best_responded_cost(p, λ1, λ2, α, μ, ν)
Returns the cost given an input p assuming that the opponent will best respond
to this p. Calculated analytically, more performant than `best_responded_cost`.
"""
function analytic_best_responded_cost(p, λ1, λ2, α, μ, ν)
    real(-(p - 1)/(μ(α) - λ1 + λ1*p) - p/(λ1*p - ν(α) + (λ2*(λ2*Complex(μ(α)*(ν(α) - λ1*p))^(1/2) - μ(α)*Complex(μ(α)*(ν(α) - λ1*p))^(1/2) - ν(α)*Complex(μ(α)*(ν(α) - λ1*p))^(1/2) - λ2*ν(α) + 2*μ(α)*ν(α) + λ1*p*Complex(μ(α)*(ν(α) - λ1*p))^(1/2) + λ1*λ2*p - 2*λ1*μ(α)*p))/(λ2*μ(α) - λ2*ν(α) + λ1*λ2*p)))
end

"""
    ∂best_responded_cost_∂p(p, λ1, λ2, α, μ, ν)
Returns the derivative of the `best_responded_cost` using AD, just as accurate
but less performant than `diff_best_responded_cost`
"""
function ∂best_responded_cost_∂p(p, λ1, λ2, α, μ, ν)
    ForwardDiff.derivative(ϕ -> best_responded_cost(ϕ, λ1, λ2, α, μ, ν), p)
end


"""
    diff_best_responded_cost(p, λ1, λ2, α, μ, ν)
Returns the derivative of the `best_responded_cost` calculated analytically.
More performant than `∂best_responded_cost_∂p`
"""
function diff_best_responded_cost(p, λ1, λ2, α, μ, ν)
    squareroot_factor = real(Complex(μ(α)*(ν(α) - λ1*p))^(1/2))
    A = (
        (λ1*(p - 1)) / (μ(α) - λ1 + λ1*p)^2 -
        (μ(α) - ν(α) + λ1*p)/
        ((squareroot_factor - ν(α) + λ1*p)*(λ2 - μ(α) - ν(α) + λ1*p)) -
        1/(μ(α) - λ1 + λ1*p)
    )  
    B = (
        p * λ1 *(
            4*real(Complex(μ(α)*(ν(α) - λ1*p))^(3/2)) -
            μ(α)^3 + 
            λ2*μ(α)^2 - 
            2*λ1^2*p^2*squareroot_factor -
            4*μ(α)^2*(ν(α) - λ1*p) +
            2*μ(α)^2*squareroot_factor -
            2*ν(α)^2*squareroot_factor -
            λ1*μ(α)*p*(ν(α) - λ1*p) +
            4*λ1*ν(α)*p*squareroot_factor +
            λ2*μ(α)*(ν(α) - λ1*p) +
            μ(α)*ν(α)*(ν(α) - λ1*p) -
            2*λ2*μ(α)*squareroot_factor
        )
    )

    C = (
        2*squareroot_factor *
        (
            squareroot_factor -
            ν(α) +
            λ1*p
        )^2 *
        (λ2 - μ(α) - ν(α) + λ1*p)^2
    )
    return real(A - B/C)
end

"""
    old_equilibrium_probability(λ1, λ2, α, μ, ν; tolerance = 10e-9, first_actor = 1, returnall = false)
Superceded by the current `equilibrium_probability` function which is much faster.
Returns the equilibrium probabilities (i.e. Nash equilibrium) for a given input.
This is calculated by finding an initial minimiser of the cost function and
recursively allowing players to best respond to the other player's actions,
stopping after a certain tolerance level is met
It has optional kwargs:
* tolerance - specifies the absolute difference required to stop the recursion
* first_actor - specifies which player should act first (this may or may not
give different results)
* returnall - defaults to false. If true return the full recursive sequence
"""
function old_equilibrium_probability(λ1, λ2, α, μ, ν; tolerance = 10e-13, first_actor = 1, returnall = false)
    if first_actor == 1
        parameters = (λ1, λ2, α, μ, ν)
        swapped_parameters = (λ2, λ1, α, μ, ν)
    else
        parameters = (λ2, λ1, α, μ, ν)
        swapped_parameters = (λ1, λ2, α, μ, ν)
    end
    
    root, ismin = classify_roots(p -> diff_best_responded_cost(p, parameters...), find_valid_roots(p -> diff_best_responded_cost(p, parameters...)))
    if sum(ismin) == 0
        val_at_zero = analytic_best_responded_cost(0,  parameters...)
        val_at_one = analytic_best_responded_cost(1,  parameters...)
        # println("it changed something")
        root = findmin((val_at_zero, val_at_one))[2]
    else
        root = first(root[ismin])
    end
    first_pstar(p) = pstar(p, parameters...)
    second_pstar(p) = pstar(p, swapped_parameters...)

    diff = 100
    response1 = 0
    ps = Float64[]
    qs = Float64[]
    if tolerance < 10e-13
        # println("tolerance of less than 10e-13 is dangerous")
        tolerance = 10e-14
    end
    while diff > tolerance
        response1 = second_pstar(root)
        response2 = first_pstar(response1)
        push!(qs, response1)
        push!(ps, response2)
        diff = abs(response2 - root)
        root = response2
    end
    if returnall == true
        return ps, qs
    else
        return root, response1
    end
end

"""
    equilibrium_probability(λ1, λ2, α, μ, ν)
Returns the equilibrium probabilities given parameters by solving the recurrence
relation p*(q*(p_n)) = p_{n+1}
"""
function equilibrium_probability(λ1, λ2, α, μ, ν)
    rec_rel(p) = pstar(qstar(p, λ1, λ2, α, μ, ν), λ1, λ2, α, μ, ν)
    p, _ = find_valid_roots(p -> rec_rel(p) - p)
    p = first(p)
    q = qstar(p, λ1, λ2, α, μ, ν)
    return (p, q)
end

function shaded_regions(mask::BitVector, iterator)
    values, indices = rle(mask)
    cum_indices = cumsum(indices)
    intervals = repeat([[0.0, 0.0]], length(indices))
    intervals[1] = iterator[[1, indices[1]]]
    for i in 2:length(cum_indices)
        intervals[i] = iterator[[cum_indices[i-1], cum_indices[i]]]
    end
    return intervals[.!values]
end

function shapes_from_intervals(intervals)
    shapes = Vector{Shape}(undef, length(intervals))
    for (i, interval) in enumerate(intervals)
        shapes[i] = Shape([(interval[1], 0), (interval[end], 0), (interval[end], 1), (interval[1], 1)])
    end
    return shapes
end

function resize_interval_shapes(
    vector_of_vectors_of_shapes::Vector{Vector{Shape}};
    n = length(vector_of_vectors_of_shapes),
    baseheight = 0)

    new_vec_vec_shapes = Vector{Vector{Shape}}(undef, length(vector_of_vectors_of_shapes))
    for i in 1:length(vector_of_vectors_of_shapes)
        new_shapes = Vector{Shape}(undef, length(vector_of_vectors_of_shapes[i]))
        for (j, shape) in enumerate(vector_of_vectors_of_shapes[i])
            new_ys = (shape.y .* (1/n)) .+ (i-1)/n .+ baseheight
            new_shapes[j] = Shape(shape.x, new_ys)
        end
        new_vec_vec_shapes[i] = new_shapes
    end
    return new_vec_vec_shapes
end

# function resize_interval_shapes(vector_of_shapes::AbstractVector{Shape})
#     n = length(vector_of_shapes)
#     new_shapes = Vector{Shape}(undef, length(vector_of_shapes))
#     for (j, shape) in enumerate(vector_of_shapes)
#         new_ys = (shape.y * (1/n)) .+ (j-1)/n
#         new_shapes[j] = Shape(shape.x, new_ys)
#     end
#     return new_shapes
# end

function calculate_assumptions(ps, qs, λ1, λ2, α, μ, ν; which = "both")
    if which == "mu" || which == "μ"
        return μ.(α) .> λ1 .* (1 .- ps)
    elseif which == "nu" || which == "ν"
        return ν.(α) .> λ1 .* ps .+ λ2 .* qs
    else
        return (μ.(α) .> λ1 .* (1 .- ps), ν.(α) .> λ1 .* ps .+ λ2 .* qs)
    end
end

function plot_assumptions!(plot, ps_list, qs, λ1, λ2, α, μ, ν;
    labels = ["μ violated", "ν violated", "a_μ violated", "a_ν violated"],
    colours = [:red, :blue, :yellow, :gray],
    height = 1/length(ps_list),
    baseheight = 0)

    assumptions = []
    for ps in ps_list
        append!(assumptions, [calculate_assumptions(ps, qs, λ1, λ2, α, μ, ν)...])
    end

    vec_vec_shapes = map(
        assumption -> shapes_from_intervals(shaded_regions(assumption, qs)),
        assumptions
    )

    new_vec_vec_shapes = resize_interval_shapes(vec_vec_shapes; n = 1/height, baseheight = baseheight)

    for (i, vec_shapes) in enumerate(new_vec_vec_shapes)
        plot!(plot, vec_shapes, colour=colours[i], label = [labels[i] fill(nothing, 1, length(vec_shapes)-1)], alpha=0.3)
    end
    return plot
end

function dp_string(number::Real, n::Integer)
    rpad(string(round(number, sigdigits = n)), n+1, "0")[1:n+1]
end

## Homo-Egualis-related functions ##
he_term(p, q, λ1, λ2, α, μ, ν, β) = β/2 * (cost(q, p, λ2, λ1, α, μ, ν) - cost(p, q, λ1, λ2, α, μ, ν))
he_cost(p, q, λ1, λ2, α, μ, ν, β) = cost(p, q, λ1, λ2, α, μ, ν) + he_term(p, q, λ1, λ2, α, μ, ν, β)

function he_pstar(q, λ1, λ2, α, μ, ν, β)
    mu = μ(α)
    nu = ν(α)
    A = (
            (mu + nu - λ1 - λ2*q)*real(Complex(-mu*(β - 2)*(2*nu - β*nu - 2*λ2*q + β*λ1*q + β*λ2*q))^(1/2)) +
            2*(
                λ1*nu +
                -2*mu*nu +
                -λ1*λ2*q +
                2*λ2*mu*q
            ) +
            β*(-λ1*nu + 2*mu*nu + λ1^2*q + λ1*λ2*q -λ1*mu*q - 2*λ2*mu*q)
        )
    B = λ1*((β-2)*(mu - nu) + q*(β*(λ1 + λ2) - 2*λ2))
    return clamp(A/B, 0, 1)
end

he_qstar(p, λ1, λ2, α, μ, ν, β) = he_pstar(p, λ2, λ1, α, μ, ν, β)

function he_equilibrium_probability(λ1, λ2, α, μ, ν, β)
    rec_rel(p) = he_pstar(he_qstar(p, λ1, λ2, α, μ, ν, β), λ1, λ2, α, μ, ν, β)
    p, _ = find_valid_roots(p -> rec_rel(p) - p)
    p = first(p)
    q = he_qstar(p, λ1, λ2, α, μ, ν, β)
    return (p, q)
end

"""
    poa(λ1, λ2, α, μ, ν)
Calculates the Price of Anarchy (ratio of worst equilibrium cost to best
possible cost), returns Inf if infeasible
"""
function poa(λ1, λ2, α, μ, ν)
    parameters = (λ1, λ2, α, μ, ν)
    he_parameters = (λ1, λ2, α, μ, ν, 1)
    tc_he = total_cost(he_equilibrium_probability(he_parameters...)..., parameters...)
    tc = total_cost(equilibrium_probability(parameters...)..., parameters...)
    ca(x, y) = calculate_assumptions(
        equilibrium_probability(x, y, α, μ, ν)...,
        x, y, α, μ, ν;
        which = "both"
    )
    return all(ca(λ1, λ2)) ? tc/tc_he : Inf
end
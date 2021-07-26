### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 3b98e7cd-b349-43c9-b45d-3b12d59bb096
begin
	using ForwardDiff
	using Roots
	using StatsBase
	using Plots
	
	
	function private_cost(p, λ1, α, μ)
	    1 / (μ(α) - λ1 + λ1 * p)
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
	    function pstar(q, λ1, λ2, α, μ, ν)
	Returns the best response to an opponents move from the perspective of the
	player responding i.e. the argument λ1  may be λ2 in practice if player 2 is
	responding to player 1. Calculated numerically, less accurate than `analytic_pstar`
	"""
	function pstar(q, λ1, λ2, α, μ, ν)
	    roots = find_valid_roots(p -> p - diff_cost(p, q, λ1, λ2, α, μ, ν))
	    if length(roots[1]) > 1
	        root_values, root_is_minimum =
	            classify_roots(p -> diff_cost(p, q, λ1, λ2, α, μ, ν), roots)
	        return (root_values[root_is_minimum])[1]
	    else
	        return (roots[1])[1]
	    end
	end
	
	
	"""
	    analytic_pstar(q, λ1, λ2, α, μ, ν)
	Returns the best response to an opponents move from the perspective of the
	player responding i.e. the argument λ1  may be λ2 in practice if player 2 is
	responding to player 1. Calculated analytically, more accurate than `pstar`
	"""
	function analytic_pstar(q, λ1, λ2, α, μ, ν)
	    p = (λ1*Complex(μ(α)*(ν(α) - λ2*q))^(1/2) - μ(α)*Complex(μ(α)*(ν(α) - λ2*q))^(1/2) - ν(α)*Complex(μ(α)*(ν(α) - λ2*q))^(1/2) - λ1*ν(α) + 2*μ(α)*ν(α) + λ2*q*Complex(μ(α)*(ν(α) - λ2*q))^(1/2) + λ1*λ2*q - 2*λ2*μ(α)*q)/(λ1*μ(α) - λ1*ν(α) + λ1*λ2*q) |> real
	    return clamp(p, 0, 1)
	end
	
	# function analytic_pstar(q, λ1, λ2, α, μ, ν)
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
	    cost(p, analytic_pstar(p, λ2, λ1, α, μ, ν), parameters...)
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
	    equilibrium_probability(λ1, λ2, α, μ, ν; tolerance = 10e-9, first_actor = 1, returnall = false)
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
	function equilibrium_probability(λ1, λ2, α, μ, ν; tolerance = 10e-9, first_actor = 1, returnall = false)
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
	    first_pstar(p) = analytic_pstar(p, parameters...)
	    second_pstar(p) = analytic_pstar(p, swapped_parameters...)
	
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
	        push!(ps, response2)
	        push!(qs, response2)
	        diff = abs(response2 - root)
	        root = response2
	    end
	    if returnall == true
	        return ps, qs
	    else
	        return root, response1
	    end
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
end

# ╔═╡ a6e91cd2-80f9-4ed4-9592-75f606001c3f
using PlutoUI

# ╔═╡ 30697798-aeac-4568-9bea-8edf7eb464bb
md"Player 2 arrival rate (λ2):  $(@bind λ2 Slider(0.00 : 0.01 : 2, default = 1, show_value = true))"

# ╔═╡ b659455c-7142-477d-b3c4-8a8e2ef56465
md"Player 1 arrival rate (λ1):  $(@bind λ1 Slider(0.00 : 0.01 : 2, default = 1, show_value = true))"

# ╔═╡ 08a1f7b9-9d03-45e9-a1dc-8972dc94602f
md"Job size (α):  $(@bind α Slider(0.00 : 0.01 : 2, default = 1, show_value = true))"

# ╔═╡ c08b6678-6b50-472a-b236-4e5c718e1aac
begin
    μ(α) = 1 / (1.5 * (α))
    ν(α) = 1 / (0.7 * (α))
    qs = 0:0.001:1
end

# ╔═╡ 7f9cc5b5-9a74-49b5-8798-dc9003f5ed7b
begin
	p, q = equilibrium_probability((λ1, λ2, α, μ, ν)...)
end

# ╔═╡ 1fbaf942-bd3f-498d-bc97-549bb70f2bce
md"
Optimal probability for Player 1:$(round(p, digits=3)),

Optimal probability for Player 2:$(round(q, digits=3))
"

# ╔═╡ 8b295584-3bfc-4b2c-bf0c-aaae22e20052
begin
    plot(qs, private_cost.(p, qs, (λ1, λ2, α, μ, ν)...), label = "player 1 private cost")
    plot!(qs, private_cost.(qs, p, (λ2, λ1, α, μ, ν)...), label = "player 2 private cost")
    plot!(qs, public_cost.(p, qs, (λ1, λ2, α, μ, ν)...), label ="shared public cost")
    hline!([0], lc = :black, label = "")
    # vline!([p], label = "p*")
    # vline!([q], label = "q*")
    ylims!(-100,100)
	plot!(legend = :outerbottom, size = (600,600))
	plot!(xlabel = "Player 2 probability (q)",
		ylabel = "Cost")
    # xlims!(0.4, 0.7)
end

# ╔═╡ 4ff08cec-ea10-45e0-8d1d-370327a16fa1
md"
### *Why it breaks*
Here we can inspect why exactly a given setup breaks. By changing `parms` we can see which of the servers becomes overloaded.
"

# ╔═╡ bf75f2d2-5f8b-42e9-a68c-4fc3cd0f2201
parms = (3, 3, 0.5, μ, ν)

# ╔═╡ 127d39cc-03ca-4bcf-b8a5-260b15a42115
md"Here are the optimal probabilities for each player"

# ╔═╡ bb9b6889-09a1-4496-bf5e-77a09c9f94bc
pr, qr = equilibrium_probability(parms...)

# ╔═╡ 5782b462-9161-4ac3-a82a-8a09e029ec14
md"and here we can compare the servers for player 1:"

# ╔═╡ 68aaf196-76f5-4348-a32e-1e85c8a3a99e
md"Total cost for Player 1: $(round(cost(pr, qr, parms...), digits=3))"

# ╔═╡ a9febb49-bf0e-48bf-a3fc-8c7cda7f14a2
round(pr*parms[1]*parms[3], digits=3)

# ╔═╡ 49a4411d-d740-4db0-8fee-e8e8ad081e76
md"Maximum possible public arrival rate: $(round(ν(parms[3]), digits=3))"

# ╔═╡ bb6b04a5-c95e-4c88-b77d-43303fc41d07
round((1-pr)*parms[1]*parms[3],digits=3)

# ╔═╡ 8925b3f7-4b74-4503-8cbc-e7e2a48c8f1e
round(μ(parms[3]),digits=3)

# ╔═╡ Cell order:
# ╟─3b98e7cd-b349-43c9-b45d-3b12d59bb096
# ╟─a6e91cd2-80f9-4ed4-9592-75f606001c3f
# ╟─30697798-aeac-4568-9bea-8edf7eb464bb
# ╟─b659455c-7142-477d-b3c4-8a8e2ef56465
# ╟─08a1f7b9-9d03-45e9-a1dc-8972dc94602f
# ╟─c08b6678-6b50-472a-b236-4e5c718e1aac
# ╟─7f9cc5b5-9a74-49b5-8798-dc9003f5ed7b
# ╟─1fbaf942-bd3f-498d-bc97-549bb70f2bce
# ╠═8b295584-3bfc-4b2c-bf0c-aaae22e20052
# ╠═4ff08cec-ea10-45e0-8d1d-370327a16fa1
# ╠═bf75f2d2-5f8b-42e9-a68c-4fc3cd0f2201
# ╟─127d39cc-03ca-4bcf-b8a5-260b15a42115
# ╠═bb9b6889-09a1-4496-bf5e-77a09c9f94bc
# ╟─5782b462-9161-4ac3-a82a-8a09e029ec14
# ╠═68aaf196-76f5-4348-a32e-1e85c8a3a99e
# ╠═a9febb49-bf0e-48bf-a3fc-8c7cda7f14a2
# ╠═49a4411d-d740-4db0-8fee-e8e8ad081e76
# ╠═bb6b04a5-c95e-4c88-b77d-43303fc41d07
# ╠═8925b3f7-4b74-4503-8cbc-e7e2a48c8f1e

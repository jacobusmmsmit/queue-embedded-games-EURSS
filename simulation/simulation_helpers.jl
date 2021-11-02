using Plots
using DataFrames
using Distributions
using ShiftedArrays
using CategoricalArrays

function simulate_queue(
    until::Real,
    interval_length::Real,
    # job_size_dist::Distribution,
    p1_arrival_dist::Distribution,
    p2_arrival_dist::Distribution,
    service_dist_pub::Distribution,
    service_dist_priv::Distribution,
    p1_prob_public::Real,
    p2_prob_public::Real,)
    now = 0.0
    p1_next_arrival = rand(p1_arrival_dist)
    p2_next_arrival = rand(p2_arrival_dist)
    p1_priv_next_departure = Inf
    p2_priv_next_departure = Inf
    public_next_departure = Inf
    events = [
        p1_next_arrival,
        p2_next_arrival,
        p1_priv_next_departure,
        p2_priv_next_departure,
        public_next_departure,
    ]
    p1_queue = Float64[]
    p2_queue = Float64[]
    public_queue = Float64[]
    public_queue_who = Int8[]
    event_data = DataFrame(time = Union{Float64,Missing}[], eventtype = String[])
    interval_data = DataFrame(interval = Int[], p1_pp = Float64[], p2_pp = Float64[])
    interval = 1
    while now < until
        push!(interval_data, [interval, p1_prob_public, p2_prob_public])
        while now < interval_length * (interval)
            nextevent = argmin(events)
    # p1_next_arrival
            if nextevent == 1
                now = p1_next_arrival
                where_to = rand(Bernoulli(p1_prob_public))
                if where_to == 1
                    if isempty(public_queue)
                        public_next_departure = p1_next_arrival + rand(service_dist_pub)
                    end
                    push!(public_queue, p1_next_arrival)
                    push!(public_queue_who, 1)
                    push!(event_data, [now "p1-public_arrival"])
                else
                    if isempty(p1_queue)
                        p1_priv_next_departure = p1_next_arrival + rand(service_dist_priv)
                    end
                    push!(p1_queue, p1_next_arrival)
                    push!(event_data, [now "p1-private_arrival"])
                end
                p1_next_arrival += rand(p1_arrival_dist)

    # p2_next_arrival
            elseif nextevent == 2
                now = p2_next_arrival
                where_to = rand(Bernoulli(p2_prob_public))
                if where_to == 1
                    if isempty(public_queue)
                        public_next_departure = p2_next_arrival + rand(service_dist_pub)
                    end
                    push!(public_queue, p2_next_arrival)
                    push!(public_queue_who, 2)
                    push!(event_data, [now "p2-public_arrival"])
                else
                    if isempty(p2_queue)
                        p2_priv_next_departure = p2_next_arrival + rand(service_dist_priv)
                    end
                    push!(p2_queue, p2_next_arrival)
                    push!(event_data, [now "p2-private_arrival"])
                end
                p2_next_arrival += rand(p2_arrival_dist)

                # p1_priv_next_departure
                elseif nextevent == 3
                    now = p1_priv_next_departure
                    push!(event_data, [now "p1-private_departure"])
                    popfirst!(p1_queue)
                    p1_priv_next_departure += isempty(p1_queue) ? Inf : rand(service_dist_priv)

                # p2_priv_next_departure
                elseif nextevent == 4
                    now = p2_priv_next_departure
                    push!(event_data, [now "p2-private_departure"])
                    popfirst!(p2_queue)
                    p2_priv_next_departure += isempty(p2_queue) ? Inf : rand(service_dist_priv)

                # public_next_departure
                else
                    now = public_next_departure
                    push!(event_data, [now "p$(popfirst!(public_queue_who))-public_departure"])
                    popfirst!(public_queue)
                    public_next_departure += isempty(public_queue) ? Inf : rand(service_dist_pub)
            end
            events = [
                p1_next_arrival,
                p2_next_arrival,
                p1_priv_next_departure,
                p2_priv_next_departure,
                public_next_departure,
            ]
        end
        interval += 1
    # choose new probabilities (i.e. strategies) here
        # p1_prob_public = p1_prob_public
        # p2_prob_public = p2_prob_public
    end
    return event_data[1:end - 1, :], interval_data
end

"""
simulate_queue

Alternative version where you input only numbers
"""
function simulate_queue(sim_length::Real,
                        int_length::Real,
                        p1_AR::Real,
                        p2_AR::Real,
                        pu_SR::Real,
                        pr_SR::Real,
                        js::Real,)
    processing_rate(base_rate, job_size) = base_rate * (1 / job_size)
    pst, qst = equilibrium_probability(
        p1_AR,
        p2_AR,
        js,
        α -> processing_rate(pr_SR, α),
        α -> processing_rate(pu_SR, α)
    )
    
    simulate_queue(
        sim_length,
        int_length,
        Exponential(1 / p1_AR),
        Exponential(1 / p2_AR),
        Exponential(processing_rate(pu_SR, js)),
        Exponential(processing_rate(pr_SR, js)),
        pst,
        qst
    )    
end

function prepare_event_data_for_plot(event_data::DataFrame, colours::AbstractVector = [:red, :orange, :blue, :green])
    @assert length(colours) == 4
    arr_dept_to_effect(x) = x == "arrival" ? 1 : -1
    event_plot_data = copy(event_data)
    event_plot_data[:, :effect] = arr_dept_to_effect.(event_plot_data[:, :event])
    sort!(event_plot_data, :time)
    event_plot_data = transform!(groupby(event_plot_data, :source),
        :effect => cumsum,
        :time => Base.Fix2(lag, -1) => :changetime,
        ungroup = false
    )
    event_plot_data = transform!(event_plot_data, :effect_cumsum => Base.Fix2(lag, -1) => :next_effect_cumsum, ungroup = true)

    colour_dict = Dict("p1-private" => colours[1], "p1-public" => colours[2], "p2-private" => colours[3], "p2-public" => colours[4])
    translate(x) = colour_dict[x]
    event_plot_data[:, :groupcolour] = translate.(event_plot_data[!, :source])
    return event_plot_data
end

function plot_event_data(plottable_event_data::DataFrame)
    p = plot()
    for g in groupby(plottable_event_data, :source)
        plot!(p, g.time, g.effect_cumsum, colour = first(g.groupcolour), label = first(g.source))
    end
    return p
end

# plot_event_data(prepare_event_data_for_plot(event_data))

function plot_path(df)
    p = plot(title = "Simulation of a queue", xlabel = "Time", ylabel = "Queue length")
    @inbounds for row in eachrow(df)
        plot!(p,
            [row.time, row.changetime],
            [row.effect_cumsum, row.effect_cumsum],
            colour = row.groupcolour,
            label = false
        )
        plot!(p,
            [row.changetime, row.changetime],
            [row.effect_cumsum, row.next_effect_cumsum],
            colour = row.groupcolour,
            label = false)
    end
    firstlinesdf = combine(groupby(df, :source)) do sdf
        sdf[argmin(skipmissing(sdf.time)), :]
    end
    @inbounds for row in eachrow(firstlinesdf)
        plot!(p,
            [0, row.time],
            [0, 0],
            colour = row.groupcolour,
            label = row.source
        )
    end

    tempdf = transform(df[:, [:time, :effect]], :effect => cumsum, :time => Base.Fix2(lag, -1) => :changetime)
    transform!(tempdf, :effect_cumsum => Base.Fix2(lag, -1) => :next_effect_cumsum)
    @inbounds for row in eachrow(tempdf)
        plot!(p,
            [row.time, row.changetime],
            [row.effect_cumsum, row.effect_cumsum],
            colour = :black,
            label = false
        )
        plot!(p,
            [row.changetime, row.changetime],
            [row.effect_cumsum, row.next_effect_cumsum],
            colour = :black,
            label = false,
            legend = :topleft)
    end
    return p
end

# p = plot_path(prepare_event_data_for_plot(event_data))

# savefig(p, "outputs/simulated_queue_plot.pdf")



function match_events(df::AbstractDataFrame)
    source = df[1, :source]
    arrivals, departures = groupby(df, :event)
    arrivals = DataFrame(arrivals)
    departures = DataFrame(departures)
    # Fill rows in departures so that it has the same length a arrivals
    for _ in 1:(nrow(arrivals) - nrow(departures))
        push!(departures, [missing source "departures"])
    end
    rename!(departures, :time => :departure_time)
    rename!(arrivals, :time => :arrival_time)
    select!(departures, Not([:source, :event]))
    # Concatenate arrivals and departures so that the n-th arrival is on the same row as the n-th departure
    hcat(departures, arrivals, makeunique=true)
end

"""
Warning: Will give missing if there are no arrivals of either private or public for a certain player in a certain time
"""
function evaluate_cost(df)
    if nrow(df) != 2
        return missing
    end
    private_wait_time = df[2, :mean_response_time]
    public_wait_time = df[1, :mean_response_time]
    player = df[1, :player]
    p = first(df[1, Regex("$(player)")])
    p * public_wait_time + (1 - p) * private_wait_time
end
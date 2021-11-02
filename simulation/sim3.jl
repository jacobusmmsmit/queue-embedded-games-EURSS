using StatsPlots

include("simulation_helpers.jl")
include("../helper_functions.jl")

begin
    begin
        total_simulation_length = 1000
        interval_length = 1000
        job_size = 1
        p1_arrival_rate = 1
        p2_arrival_rate = 1
        private_service_rate = 1 / 1.5
        public_service_rate = 1 / 0.7
        breakin = 3
        pst, qst = equilibrium_probability(
        p1_arrival_rate,
        p2_arrival_rate,
        job_size,
        α -> 1 / (1.5 * (α)),
        α -> 1 / (0.7 * (α)))
    end

event_data, interval_data = simulate_queue(
    total_simulation_length,
    interval_length,
    Exponential(1 / p1_arrival_rate),
    Exponential(1 / p2_arrival_rate),
    Exponential(job_size / public_service_rate),
    Exponential(job_size / private_service_rate),
    pst,
    qst,
)

transform!(event_data, :eventtype => ByRow(x -> split(x, "_")) => [:source, :event])
select!(event_data, Not(:eventtype))

event_data_wide = vcat([match_events(gdf) for gdf in groupby(event_data, :source)]...)
select!(event_data_wide, [:source, :arrival_time, :departure_time])
transform!(event_data_wide,
    [:arrival_time, :departure_time] => ((x, y) -> (y - x)) => :response_time,
    :arrival_time => (x -> cut(x, 0:interval_length:total_simulation_length + interval_length, extend=true)) => :interval)
performance_data = combine(groupby(event_data_wide, [:interval, :source]), :response_time => (x -> mean(skipmissing(x))) => :mean_response_time)

transform!(performance_data, :source => ByRow(x -> split(x, "-")) => [:player, :dummy])
select!(performance_data, Not(:dummy))
groupby(performance_data, [:player, :interval])
interval_data.cut = cut((interval_data.interval .- 1) .* interval_length, 0:interval_length:total_simulation_length + interval_length)
performance_data = leftjoin(performance_data, select(interval_data, Not(:interval)), on = [:interval => :cut])

cost_data = combine(groupby(performance_data, [:player, :interval])) do sdf
    evaluate_cost(sdf)
end
        rename!(cost_data, :x1 => :cost)

begin
    p2 = plot(title = "Cost by Communication Interval", xlabel = "Interval", ylabel = "Cost")
    cost_data = leftjoin(sort!(cost_data, :player), interval_data[:, [:interval, :cut]], on = (:interval => :cut), makeunique=true)
    for g in groupby(cost_data, :player)
        plot!(p2, g.interval_1, g.cost, label = first(g.player), lw = 2)
    end
    p2
end
p2 
end

    println(median(skipmissing(event_data_wide.response_time)))
histogram(log10.(collect(skipmissing(event_data_wide.response_time))),
    nbins = 50,
    label = false,
    xlabel = "log10(Response Time)",
    ylabel = "Frequency",
    title = "Histogram of Response Times")


ARs = 0.1:0.1:2
total_volume = 2
JSs = 1 / 2 * total_volume ./ ARs
n_reps = 50
dot_every = Int(floor(n_reps / 5))

begin
mean_resval = Float64[]
var_resval = Float64[]

for (AR, JS) in zip(ARs, JSs)
    mean_temp = Float64[]
    var_temp = Float64[]
    private_service_rate = 1 / (1.5 * JS)
    public_service_rate = 1 / (0.7 * JS)

    print("AR = ", AR, " JS = ", JS, " ")
    for i in 1:n_reps
        if mod(i, dot_every) == 0
            print(".")
        end
        event_data, interval_data = simulate_queue(
            total_simulation_length,
            total_simulation_length,
            AR,
            AR,
            1 / 0.7,
            1 / 1.5,
            JS
        )
        transform!(event_data, :eventtype => ByRow(x -> split(x, "_")) => [:source, :event])
        select!(event_data, Not(:eventtype))
        event_data_wide = vcat([match_events(gdf) for gdf in groupby(event_data, :source)]...)
        select!(event_data_wide, [:source, :arrival_time, :departure_time])
            transform!(event_data_wide,
                    [:arrival_time, :departure_time] => ((x, y) -> (y - x)) => :response_time)
        mean_and_var = combine(event_data_wide,
        :response_time => (
            x -> (mean_RT = mean(skipmissing(x)), var_RT = var(skipmissing(x)))
            ) => AsTable
        )

        append!(mean_temp, mean_and_var[1, 1])
        append!(var_temp, mean_and_var[1, 2])
    end
    println()
    append!(mean_resval, mean(mean_temp))
    append!(var_resval, mean(var_temp))
end
end

plot(ARs, JSs, xlabel = "Arrival Rate", ylabel = "Job Size")
mean_plot = plot(JSs, mean_resval, xlabel = "Job Size", yaxis = :log10, ylabel = "Mean Response Time", legend = :none)
sd_plot = plot(JSs, sqrt.(var_resval), xlabel = "Job Size", yaxis = :log10, ylabel = "SD Response Time", legend = :none)

begin
    total_simulation_length = 10000
    interval_length = total_simulation_length / 10
    job_size = 0.5
    p1_arrival_rate = 2
    p2_arrival_rate = 2
    private_service_rate = 1 / (1.5 * job_size)
    public_service_rate = 1 / (0.7 * job_size)
    breakin = 3
    pst, qst = equilibrium_probability(
        p1_arrival_rate,
        p2_arrival_rate,
        job_size,
        α -> 1 / (1.5 * (α)),
        α -> 1 / (0.7 * (α)))
end

# to ensure correctness: calculate_distributions

event_data, interval_data = simulate_queue(
    1000,
    1000,
    1,
    1,
    1 / 0.7,
    1 / 1.5,
    1
)

function evaluate_response_time(event_data)
    transform!(event_data, :eventtype => ByRow(x -> split(x, "_")) => [:source, :event])
    select!(event_data, Not(:eventtype))
    event_data_wide = vcat([match_events(gdf) for gdf in groupby(event_data, :source)]...)
    select!(event_data_wide, [:source, :arrival_time, :departure_time])
    transform!(event_data_wide,
                [:arrival_time, :departure_time] => ((x, y) -> (y - x)) => :response_time)
    df = combine(event_data_wide,
        :response_time => (
            x -> (mean_RT = mean(skipmissing(x)), sd_RT = sqrt(var(skipmissing(x))))
        ) => AsTable
    )
    return df[1, 1], df[1, 2]
end

evaluate_response_time(event_data)

prepared_data = prepare_event_data_for_plot(event_data)
@df prepared_data plot(
    :time,
    :effect_cumsum,
    group = :groupcolour,
)
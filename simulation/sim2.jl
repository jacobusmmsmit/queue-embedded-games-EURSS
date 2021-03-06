using Random
using StatsBase
using Distributions
using StatsPlots
using DataFrames
using ShiftedArrays
using Measures

Random.seed!(123)

begin
    arrival_rates = (0.1, 0.1)
    job_size = 10
    private_service_rate(job_size) = 1 / (1.5 * job_size)
    public_service_rate(job_size) = 1 / (0.7 * job_size)
    p1_prob_public = 0.576
    p2_prob_public = 0.576
    segment_time = 10000
    total_time = 10000
end

function run_simulation(arrival_rates, job_size, private_service_rate, public_service_rate, p1_prob_public, p2_prob_public; total_time = 1000, return_ft = false)
    function generate_interarrival_times(d::Distribution, time)
        interarrival_times = [rand(d)]
            # initial_n = time/mean(d)
            # append!(interarrival_times, rand(d, initial_n))
        while sum(interarrival_times) <= time
            append!(interarrival_times, rand(d))
        end
        return interarrival_times
    end

    job_interarrival_times = generate_interarrival_times(Exponential(1 / sum(arrival_rates)), total_time)
    total_jobs = length(job_interarrival_times)
    job_which_player = 1 .+ rand(Bernoulli(arrival_rates[1] / sum(arrival_rates)), total_jobs)
    df = DataFrame(order = 1:total_jobs, interarrival_times = job_interarrival_times, player = job_which_player)
    df.public = @. ifelse(df.player == 1, rand(Bernoulli(p1_prob_public)), rand(Bernoulli(p2_prob_public)))
    df.service_duration = @. ifelse(
        df.public,
        rand(Exponential(1 / public_service_rate(job_size))),
        rand(Exponential(1 / private_service_rate(job_size)))
    )
    df.arrival_time = cumsum(df.interarrival_times)
    df.finished_time = df.arrival_time .+ df.service_duration

    # df_public = filter(a -> a.public, df)
    df_public = subset(df, :public, view = true)
    df_private = subset(df, :public => ByRow(!=(true)), view = true)
    df_private1 = subset(df_private,  :player => ByRow(==(1)), view = true)
    df_private2 = subset(df_private,  :player => ByRow(==(2)), view = true)

    function correct_finished_time!(df)
        df_return = view(df, :, :)
        @inbounds for i = 2:size(df_return, 1)
            if lag(df_return.finished_time)[i] > df_return.arrival_time[i]
                df_return[i, :finished_time] = lag(df_return.finished_time)[i] + df_return.service_duration[i]
            end
        end
        return df_return
    end

    df_public = correct_finished_time!(df_public)
    df_private1 = correct_finished_time!(df_private1)
    df_private2 = correct_finished_time!(df_private2)

    function make_events(arrivals_df)
        df = arrivals_df
        function f1(df, column, event_result, length_df = size(df)[1])
            times_df = select(df, [:player, :public, column])
            rename!(times_df, column => :time)
            times_df[!, :event] = repeat([event_result], length_df)
            return times_df
        end

        arrivaltimes_df = f1(arrivals_df, :arrival_time, 1)
        servicetimes_df = f1(arrivals_df, :finished_time, -1)
        
        return vcat(arrivaltimes_df, servicetimes_df)
    end

    if return_ft
        return select(df, [:player, :public, :arrival_time, :service_duration, :finished_time])
    end

    df_events = sort(make_events(df), :time)
    df_events_simple = select(df_events, [:time, :event])
    transform!(df_events_simple,
        :event => cumsum => :cumevent,
        :time => Base.Fix2(lag, -1) => :endtime,
    )
    df_events_simple[:, :group] .= -1
    transform!(df_events_simple,
        :cumevent => Base.Fix2(lag, -1) => :nextevent
    )
    transform!(df_events, [:player, :public] => ((x, y) -> ifelse.(y, 0, x)) => :group)
    gby = groupby(df_events, :group)
    transform!(gby,
        :event => cumsum => :cumevent,
        :time => Base.Fix2(lag, -1) => :endtime,
    )
    transform!(gby, :cumevent => Base.Fix2(lag, -1) => :nextevent)
    return(df_events_simple, gby)
end

df_events_simple, gby = run_simulation(arrival_rates, job_size, private_service_rate, public_service_rate, p1_prob_public, p2_prob_public)
gby2 = gby
gby = DataFrame(gby)
dfs = vcat(gby, df_events_simple, cols=:union)

colour_dict = Dict(-1 => :black, 0 => :red, 1 => :blue, 2 => :green)
label_dict = Dict(-1 => "Total", 0 => "Public", 1 => "Player 1 Private", 2 => "Player 2 Private")
translate(x) = colour_dict[x]
dfs[:, :groupcolour] = translate.(dfs[!, :group])

function plot_path(df, end_time)
    p = plot(title = "Simulation of a queue", xlabel = "Time", ylabel = "Queue length")
    @inbounds for i in 1:(size(df)[1] - 1)
        plot!(p,
            [df[i, :time], df[i, :endtime]],
            [df[i, :cumevent], df[i, :cumevent]],
            colour = df[i, :groupcolour],
            label = false
        )
        plot!(p,
            [df[i, :endtime], df[i, :endtime]],
            [df[i, :cumevent], df[i, :nextevent]],
            colour = df[i, :groupcolour],
            label = false)
    end
    tempdf = combine(groupby(dfs, :group)) do sdf; sdf[1, :] end
    tempdf[:, :endtime] = tempdf[:, :time]
    tempdf[:, :time] .= 0
    @inbounds for i in 1:(size(tempdf)[1])
        plot!(p,
            [tempdf[i, :time], tempdf[i, :endtime]],
            [0, 0],
            colour = tempdf[i, :groupcolour],
            label = label_dict[tempdf[i, :group]]
        )
        plot!(p,
            [tempdf[i, :endtime], tempdf[i, :endtime]],
            [0, 1],
            colour = tempdf[i, :groupcolour],
            label = false)
    end

    tempdf2 = combine(groupby(dfs, :group)) do sdf; sdf[end, :] end
    final_time = dfs[argmax(skipmissing(dfs[:, :endtime])), :endtime]
    tempdf2[:, :endtime] .= final_time
    @inbounds for i in 1:(size(tempdf2)[1])
        plot!(p,
            [tempdf2[i, :time], tempdf2[i, :endtime]],
            [0, 0],
            colour = tempdf2[i, :groupcolour],
            label = false
        )
    end
    xlims!(p, (0, end_time))
        return p
end

# p2 = plot_path(dfs, total_time)

@df dfs plot(
    :time,
    :cumevent,
    group = :groupcolour
)

begin
JSs = 0.1:0.1:1
AR_JS_ratio = 1
ARs = AR_JS_ratio ./ JSs 

# res = DataFrame(player = Int64[], public = Bool[], JS = Float64[], mean_response_time = Float64[])
resval_mean = Float64[]
resval_var = Float64[]
for (AR, JS) in zip(ARs, JSs)
    temp_mean = Float64[]
    temp_var = Float64[]
    print("AR: ", round(AR, digits = 4), ", JS: ", round(JS, digits = 4))
    num_repetitions = JS * 10
    for i in 1:num_repetitions
        if mod(i, floor(num_repetitions / 5)) == 0 print(".") end
        _, gby = run_simulation(
            AR,
            JS,
            private_service_rate,
            public_service_rate,
            p1_prob_public,
            p2_prob_public,
            total_time = 10000
        )
        gby = DataFrame(gby)
        sort!(gby, [:player, :public, :time])

        arrivals = subset(gby, :event => ByRow(isequal(1)))
        services = subset(gby, :event => ByRow(isequal(-1))) |>
            x -> select(x, [:time]) |>
            x -> rename(x, :time => :dep_time)

        gby = hcat(arrivals, services) |>
            x -> select(x, [:player, :public, :time, :dep_time]) |>
            x -> groupby(x, [:player, :public]) |>
            x -> transform(x, [:time, :dep_time] => ((x, y) -> (y - x)) => :response_time) |>
            x -> combine(x, :response_time => (x -> (mean_RT = mean(skipmissing(x)), sd_RT = var(skipmissing(x)))) => AsTable)
        append!(temp_mean, gby[1, 1])
        append!(temp_var, gby[1, 2])
    end
    println()
    append!(resval_mean, mean(temp_mean))
    append!(resval_var, mean(temp_var))
end
end
### Mean and Variance of Response Time Plots ###
begin
include("../helper_functions.jl")
function variance_of_cost(p, q, ??1, ??2, ??, ??, ??)
    (p * public_cost(p, q, ??1, ??2, ??, ??, ??))^2 + ((1 - p) * private_cost(p, ??1, ??, ??))^2
end

Hyperexponential(parameters, weights) = MixtureModel(Exponential, parameters, weights)
RTDist(p, q, ??1, ??2, ??, ??, ??) = Hyperexponential([(1/(??(??) - ??1*p - ??2*q)), (1/(??(??) - ??1*(1-p)))], [p, 1 - p])
RT = RTDist(0.57, 0.57, 1, 1, 1, private_service_rate, public_service_rate)
RTs = RTDist.(0.57, 0.57, ARs, ARs, JSs, private_service_rate, public_service_rate)
plot_mean = plot(JSs, resval_mean, xlabel = "Job Size", ylabel = "Mean Response Time", label = "Empirical Mean")
plot!(JSs,
    mean.(RTs),
    title = "Mean Response Time by Job Size",
    xlabel = "Job Size",
    ylabel = "Mean Response Time",
    label = "Theoretical Mean",
    size = [400, 300],
    legend = :topleft)

plot_variance = plot(JSs, resval_var, xlabel = "Job Size", ylabel = "Variance of Response Time", label = "Empirical Variance")
plot!(JSs,
    var.(RTs),
    title = "Variance by Job Size",
    xlabel = "Job Size",
    ylabel = "Variance of Response Time",
    label = "Theoretical Variance",
    size = [400, 300],
    legend = :topleft)

mean_var_plot = plot(plot_mean, plot_variance, size = [800, 300], margin = 12pt)
end

savefig("outputs/simulated_variance.pdf")

plot(JSs, resval_mean ./ cost.(0.57, 0.57, ARs, ARs, JSs, private_service_rate, public_service_rate))
plot(JSs, resval_var ./ variance_of_cost.(0.57, 0.57, ARs, ARs, JSs, private_service_rate, public_service_rate))

cost(0.57, 0.57, 1, 1, 1, private_service_rate, public_service_rate)
cost(0.57, 0.57, 2, 2, 0.5, private_service_rate, public_service_rate)
# sort!(gby, [:player, :public, :time])

# arrivals = subset(gby, :event => ByRow(isequal(1)))
# services = subset(gby, :event => ByRow(isequal(-1))) |>
#     x -> select(x, [:time]) |>
#     x -> rename(x, :time => :dep_time)

# hcat(arrivals, services) |>
#     x -> select(x, [:player, :public, :time, :dep_time]) |>
#     x -> groupby(x, [:player, :public]) |>
#     x -> combine(x, [:time, :dep_time] => ((x, y) -> (mean(y - x))) => :mean_response_time) |>
#     x -> combine(x, :mean_response_time => mean)[1, 1]

C = 1
??s = 0.1:0.01:5
plot(??s, ForwardDiff.derivative.(?? -> cost(0.57, 0.57, C/??, C/??, ??, private_service_rate, public_service_rate), ??s))

sim_df = run_simulation(
            1,
            1,
            private_service_rate,
            public_service_rate,
            p1_prob_public,
            p2_prob_public;
            total_time = 100000,
            return_ft = true
) |>
    x -> transform(x, [:finished_time, :arrival_time] => (-) => "response_time")


# function generate_response_times(N, p, q, ??1, ??2, ??, ??, ??)
#     pub = rand(Exponential(1/(??(??) - ??1*p - ??2*q)), N)
#     pr = rand(Exponential(1/(??(??) - ??1*(1-p))), N)
#     @. (p * pub) + ((1-p) * pr)
# end


# histogram(generate_response_times(100000, 0.57, 0.57, 1, 1, 1, private_service_rate, public_service_rate))
# cost(0.57, 0.57, 1, 1, 1, private_service_rate, public_service_rate)


density_plot = density(log.(rand(RT, 1000000)), lw = 4, label = "Theoretical Response Distribution")
density!(log.(sim_df.response_time), lw = 4, label = "Empirical Response Distribution")
plot!(xlims = [-6, 4], legend = :topleft, size = (500, 300),
    xlabel = "log(x)",
    ylabel = "f(x)",
    title = "Comparison of Densities")

# savefig("outputs/densities.pdf")

mean(RT)
mean(sim_df.response_time)
var(RT)
var(sim_df.response_time)




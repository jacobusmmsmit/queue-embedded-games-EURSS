using Random
using StatsBase
using Distributions
using StatsPlots
using DataFrames
using ShiftedArrays

Random.seed!(123)

arrival_rates = (1, 1)
# job_size_dist = DiscreteNonParametric([1], [1])
job_size = 1
private_service_rate(job_size) = 1 / (1.5 * job_size)
public_service_rate(job_size) = 1 / (0.7 * job_size)
p1_prob_public = 0.3
p2_prob_public = 0.3


p1_arrival_dist = Exponential(1 / arrival_rates[1])
p2_arrival_dist = Exponential(1 / arrival_rates[2])

private_service_dist = Exponential(1 / private_service_rate(job_size))
public_service_dist = Exponential(1 / public_service_rate(job_size))

p1_njobs = 5
p2_njobs = 5
total_njobs = p1_njobs + p2_njobs
p1_interarrival_times = rand(p1_arrival_dist, p1_njobs)
p2_interarrival_times = rand(p2_arrival_dist, p2_njobs)

p1_job_destinations = rand(Bernoulli(p1_prob_public), p1_njobs)
p2_job_destinations = rand(Bernoulli(p2_prob_public), p2_njobs)


p1_job_order = sort(DataFrame(dests = p1_job_destinations, order = 1:p1_njobs), :dests)
p2_job_order = sort(DataFrame(dests = p2_job_destinations, order = 1:p2_njobs), :dests)

function event_df(arrival_times, service_distribution, job_destinations, player, public, job_order)
    job_destinations = public ? job_destinations : map(!, job_destinations)
    filtered_job_order = filter(r -> r.dests == public, job_order)
    DataFrame(
        arrival_time = arrival_times[job_destinations],
        service_time = rand(service_distribution,
        sum(job_destinations)),
        player = player,
        public = public,
        order = filtered_job_order.order)
end

p1_public_arrivals_df = event_df(p1_interarrival_times, public_service_dist, p1_job_destinations, 1, true, p1_job_order)
p1_private_arrivals_df = event_df(p1_interarrival_times, private_service_dist, p1_job_destinations, 1, false, p1_job_order)
p2_public_arrivals_df = event_df(p2_interarrival_times, public_service_dist, p2_job_destinations, 2, true, p2_job_order)
p2_private_arrivals_df = event_df(p2_interarrival_times, private_service_dist, p2_job_destinations, 2, false, p2_job_order)

function calculate_events(public_arrivals, private_arrivals)
    function service_finishes(arrivals_df)
        sorted = sort(arrivals_df, :order)
        sorted[!, :cum_arrival_time] = cumsum(sorted.arrival_time)
        sorted[!, :service_finishes] = sorted[!, :cum_arrival_time] .+ sorted[!, :service_time]
        return sorted
    end

    arrivals_df = service_finishes(vcat(public_arrivals, private_arrivals))

    function f1(df, column, event_result, length_df = size(df)[1])
        times_df = select(df, [:player, :public, column])
        rename!(times_df, column => :time)
        times_df[!, :event] = repeat([event_result], length_df)
        return times_df
    end

    arrivaltimes_df = f1(arrivals_df, :cum_arrival_time, 1)
    servicetimes_df = f1(arrivals_df, :service_finishes, -1)
    
    return vcat(arrivaltimes_df, servicetimes_df)
end

p1_eventtimes_df = calculate_events(p1_public_arrivals_df, p1_private_arrivals_df)
p2_eventtimes_df = calculate_events(p2_public_arrivals_df, p2_private_arrivals_df)
eventtimes_df = sort(vcat(p1_eventtimes_df, p2_eventtimes_df), :time)
eventtimes_df[!, :public] = ifelse.(eventtimes_df[!, :public] .== true, "public", "private")
ped = select(deepcopy(eventtimes_df), [:time, :event])
ped[:, :event] = cumsum(ped[:, :event])
sort!(push!(ped, (0, 0)), :time)
ped[:, :endtime] = lag(ped[:, :time], -1)
ped[:, :nextevent] = lag(ped[:, :event], -1)

p = plot()
@inbounds for i = 1:(size(ped, 1) - 1)
    # plot!(p, [ped[i, :time]], [ped[i, :event]], st = :scatter, colour = :black, label = false)
    plot!(p, [ped[i, :time], ped[i, :endtime]], [ped[i, :event], ped[i, :event]], colour = :black, label = false)
    plot!(p, [ped[i, :endtime], ped[i, :endtime]], [ped[i, :event], ped[i, :nextevent]], colour = :black, label = false)
end
# plot!(p, [ped[size(ped, 1), :time]], [ped[size(ped, 1), :event]], st = :scatter, colour = :black, label = false)

ped
p

# public_njobs = sum(p1_job_destinations + p2_job_destinations)
# private_njobs = total_njobs - public_njobs
# public_service_times = rand(public_service_dist, public_njobs)
# private_service_times = rand(private_service_dist, private_njobs)

# public_service_times
# p1_interarrival_times[p1_job_destinations]
# p2_interarrival_times[p2_job_destinations]

# p1_arrivals_df = DataFrame(times = p1_interarrival_times, cumtimes = cumsum(p1_interarrival_times), source = 1, public_dest = p1_job_destinations)
# p2_arrivals_df = DataFrame(times = p2_interarrival_times, cumtimes = cumsum(p2_interarrival_times), source = 2, public_dest = p2_job_destinations)
# arrivals_df = vcat(p1_arrivals_df, p2_arrivals_df)

# public_services_df = DataFrame(service_time = public_service_times, public = true)
# private_services_df = DataFrame(service_time = private_service_times, public = false)
# services_df = vcat(public_services_df, private_services_df)
# hcat(arrivals_df, services_df)



# @df arrivals_df scatter(
#     :times,
#     1:total_njobs,
#     color = :source
# )

# plot(sort(vcat(cumsum(p1_interarrival_times), cumsum(p2_interarrival_times))), 1:20, seriestype = :scatter)
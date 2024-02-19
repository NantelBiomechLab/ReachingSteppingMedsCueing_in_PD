"Filter a marker array"
function condition_marker(sig)
    df = let
        n = 4
        bw = Butterworth(n)
        lpf = Lowpass(12.; fs=100)
        digitalfilter(lpf, bw)
    end

    return filtfilt(df, sig)
end

find_cycles(trial, src=getsource(trial, "main"); filter5=true) = find_cycles(trial, src, Val(Symbol(conditions(trial)[:task])); filter5)

struct CycleInterval{T<:AbstractRange}
    rising::T
    falling::T

    function CycleInterval(rising::T, falling::T) where {T<:AbstractRange}
        step(rising) > 0 || throw(DomainError(step(rising),
            "the `step` of `rising` must be positive"))
        step(falling) > 0 || throw(DomainError(step(falling),
            "the `step` of `falling` must be positive"))
        step(rising) == step(falling) || throw(DomainError((step(rising), step(falling)),
            "the `step` of both `rising` and `falling` ranges must be the same"))
        last(rising) == first(falling) || throw(DomainError((last(rising), first(falling)),
            "the `falling` range must overlap the `rising` range by one element"))
        return new{T}(rising, falling)
    end
end

CycleInterval(lo, mid, hi) = CycleInterval(lo:mid, mid:hi)

function Base.show(io::IO, cycle::CycleInterval)
    print(io, "CycleInterval(", first(cycle), ",", midpoint(cycle), ",", last(cycle), ")")
end

Base.isless(A::CycleInterval, B::CycleInterval) = isless(A.rising, B.rising)
Base.first(intvl::CycleInterval) = first(intvl.rising)
midpoint(intvl::CycleInterval) = last(intvl.rising)
Base.last(intvl::CycleInterval) = last(intvl.falling)
Base.collect(intvl::CycleInterval) = [intvl.rising; intvl.falling]
Base.length(intvl::CycleInterval) = length(intvl.rising) + length(intvl.falling) - 1

extents(cycle::CycleInterval) = first(cycle.rising):last(cycle.falling)


"Return an array of ranges for which y is `missing`"
function missingintervals(y)
    ints = UnitRange[]
    i = firstindex(y)
    while true
        if ismissing(y[i])
            iend = findnext(!ismissing, y, i+1) # find next existing value
            if isnothing(iend) # reached end of `y` and all are missing
                push!(ints, i:lastindex(y))
                break
            end

            push!(ints, i:(iend-1))
            inext = findnext(ismissing, y, iend+1) # find next missing value
            isnothing(inext) && break # reached end of `y` and no more missing

            i = inext
        else
            inext = findnext(ismissing, y, i+1)
            isnothing(inext) && break # reached end of `y` and no more missing

            i = inext
        end
    end

    return ints
end

"Return a histogram of the `missing`ness of `signal` among the ensemble of `cycles`"
function cycle_missing_histogram(cycles, signal; normalized=true)
    ints = missingintervals(signal)

    missed = zeros(Int,100)

    for rg in timeintervals_toindices(cycles, ints)
        missed[rg] .+= 1
    end

    if normalized
        missed ./= (length(cycles)-1)
    end
    return missed
end

"""
    contiguousranges(cycles::Vector{CycleInterval{T}}) where T -> ranges::Vector{UnitRange}}

Return an array of ranges of the contiguous `cycles`; e.g. where `ranges[i]` is a range
covering one or more cycles in `cycles` that are contiguous.
"""
function contiguousranges(cycles::Vector{CycleInterval{T}}) where T
    TI = eltype(T)
    rgs = UnitRange{TI}[]
    f = first(first(cycles))
    for i in eachindex(cycles)[1:end-1]
        if isdisjoint(extents(cycles[i]), extents(cycles[i+1]))
            push!(rgs, f:last(cycles[i]))
            f = first(cycles[i+1])
        end
    end
    push!(rgs, f:last(last(cycles)))

    return rgs
end

# Added to base Julie in PR #46356; first available in release v1.10
@static if VERSION < v"1.10"
    "Calculate if ranges `a` and `b` are disjoint"
    function Base.isdisjoint(a::AbstractRange{T}, b::AbstractRange{T}) where T
        (isempty(a) || isempty(b)) && return true
        fa, la = extrema(a)
        fb, lb = extrema(b)
        if (la < fb) | (lb < fa)
            return true
        else
            return _overlapping_range_isdisjoint(a, b)
        end
    end

    _overlapping_range_isdisjoint(a::AbstractRange{T}, b::AbstractRange{T}) where T = invoke(isdisjoint, Tuple{Any,Any}, a, b)

    function _overlapping_range_isdisjoint(a::AbstractRange{T}, b::AbstractRange{T}) where T<:Integer
        if abs(step(a)) == abs(step(b))
            return mod(minimum(a), step(a)) != mod(minimum(b), step(a))
        else
            return invoke(isdisjoint, Tuple{Any,Any}, a, b)
        end
    end
end

"Return an array of ranges where each range is `last(ranges[i]):first(ranges[i+1])`"
function invert_array_of_ranges(ranges)
    T = eltype(first(ranges))
    inv_rgs = UnitRange{T}[]
    for i in eachindex(ranges)[1:end-1]
        push!(inv_rgs, last(ranges[i]):first(ranges[i+1]))
    end

    return inv_rgs
end

"Calculate PCI (Plotnik et al. 2007) for a cycle"
function cycle_pci(cycles::Vector{CycleInterval{T}}) where T
    φ = 360 .* length.(getfield.(cycles, :rising))./length.(cycles)
    φ_ABS = mean(abs, φ .- 180)
    Pφ_ABS = (φ_ABS/180)*100
    φ_CV = variation(φ)*100

    return φ_CV + Pφ_ABS
end

function cycles_analysis(trial; analyze_cal_trial=false)
    seg = Segment(trial, "main")
    requiresource!(trial, "cal_rrt_freq" => Source{C3DFile})

    sr = SegmentResult(seg)
    res = results(sr)

    FS = 100

    # Only used in (unpublished) exploratory analysis
    if analyze_cal_trial
        cal_source = conditions(trial)[:task] == "rrt" ? "cal_rrt_freq" : "cal_sip_freq"
        cal_right_strides, cal_left_strides, cal_fog = find_cycles(trial, getsource(trial, cal_source); filter5=false)
        res["cal_fog"] = cal_fog
        target_period = median(length.(unique!([
            getfield.(cal_right_strides, :rising);
            getfield.(cal_right_strides, :falling);
            getfield.(cal_left_strides, :rising);
            getfield.(cal_left_strides, :falling);
        ])))/FS
        res["calc_target_period"] = target_period
        res["calc_target_bpm"] = 60/target_period
    end

    right_strides, left_strides, fog = find_cycles(trial, source(seg))
    all_strides = sort!([right_strides; left_strides])
    res["numfog"] = length(invert_array_of_ranges(contiguousranges(all_strides))) +
        ifelse(last(last(all_strides)) < (numpointframes(readsource(trial, "main")) - 500),
            1, 0)

    left_steps = unique!([ getfield.(right_strides, :rising); getfield.(left_strides, :falling) ])
    right_steps = unique!([ getfield.(left_strides, :rising); getfield.(right_strides, :falling) ])
    all_steps = sort!([left_steps; right_steps])
    res["numsteps"] = length(all_steps)

    # All PD group have identified ma_side, default of "right" is specified for controls,
    # but is not used in any analysis of controls, therefore is meaningless (just makes
    # the `get` function call work)
    ma_side = get(conditions(trial), :ma_side, "right")

    if ma_side == "right"
        ma = "right"
        la = "left"
    else
        ma = "left"
        la = "right"
    end

    res["left_step_time"] = mean(length, left_steps)/FS
    res["right_step_time"] = mean(length, right_steps)/FS
    res["ma_step_time"] = res["$(ma)_step_time"]*1000 # ms
    res["la_step_time"] = res["$(la)_step_time"]*1000 # ms
    res["step_time_asym"] = res["ma_step_time"] - res["la_step_time"]

    res["left_step_time_cov"] = variation(length.(left_steps)./FS, res["left_step_time"] )*100 # %
    res["right_step_time_cov"] = variation(length.(right_steps)./FS, res["right_step_time"] )*100 # %
    res["ma_step_time_cov"] = res["$(ma)_step_time_cov"]
    res["la_step_time_cov"] = res["$(la)_step_time_cov"]
    res["step_time_cov_asym"] = res["ma_step_time_cov"] - res["la_step_time_cov"]

    μ_IRI = mean(length, all_steps)/FS
    Gi_1 = mean(x -> (length(x.falling)/FS - μ_IRI)*(length(x.rising)/FS - μ_IRI), all_strides)
    Gi_0 = var(length.(all_steps)./FS; corrected=false, mean=μ_IRI)
    res["lag1_corr"] = Gi_1/Gi_0
    res["clock_var"] = Gi_0 + 2Gi_1 # s^2
    res["motor_var"] = -Gi_1 # s^2

    # res["left_step_time_err"] = (target_period - res["left_step_time"])/target_period*100 # %
    # res["right_step_time_err"] = (target_period - res["right_step_time"])/target_period*100 # %
    # res["step_time_err_asym"] = res["left_step_time_err"] - res["right_step_time_err"]

    res["pci"] = cycle_pci([left_strides; right_strides])

    res["act_bpm"] = 60/μ_IRI

    return sr
end

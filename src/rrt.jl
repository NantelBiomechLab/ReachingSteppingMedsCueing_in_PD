function condition_touch(sig)
    df = let
        n = 4
        bw = Butterworth(n)
        lpf = Lowpass(499.; fs=1000)
        digitalfilter(lpf, bw)
    end

    return filtfilt(df, sig)
end

# Memory efficient "infinite" array
struct LazyFill{T}
    lz::T
end

# Act like an infinite array
Base.getindex(l::LazyFill{T}, inds...) where T = l.lz

function find_cycles(trial, src, ::Val{:rrt}; filter5=true)
    rsrc = readsource(src; strip_prefixes=true)

    # Prioritize touch events determined by the capacitive touch target
    touchsig = condition_touch(rsrc.analog["Electric Current.DTRG"])
    if maximum(touchsig) > 2 # Target functioning correctly
        touchsig_vel = centraldiff(touchsig;
            dt=inv(rsrc.groups[:ANALOG][Float32, :RATE]))
        pks, _ = findmaxima(touchsig_vel)
        _, proms = peakproms(pks, touchsig_vel)
        peakproms!(pks, touchsig_vel; minprom=maximum(proms)/2)
        pks = round.(Int, pks ./ 10)
    else
        @debug "Touch target not connected" trial
        # Use finger markers as initial/base event source
        rmkr = @view rsrc.point["RFIN"][:,2]
        lmkr = @view rsrc.point["LFIN"][:,2]
        ridxs, = findmaxima(rmkr, 5; strict=false);
        peakproms!(ridxs, rmkr; strict=false, minprom=100)
        lidxs, = findmaxima(lmkr, 5; strict=false);
        peakproms!(lidxs, lmkr; strict=false, minprom=100)
        pks = sort!([ridxs; lidxs])
    end

    # Determine touch side (right or left)
    side = Vector{Char}(undef, length(pks))
    for (i, touch) in enumerate(pks)
        flag = 0
        for (R, L) in zip(["RIDX", "RFIN", "RWRM", "RWRL"], ["LIDX", "LFIN", "LWRM", "LWRL"])
            Rval = get(rsrc.point, R, LazyFill(missing))[max(touch, 1), 2]
            Lval = get(rsrc.point, L, LazyFill(missing))[max(touch, 1), 2]
            R_L = Rval - Lval

            # Assume that the leading hand is the reaching hand, if there is an AP
            # separation of at least 100mm or if the leading hand has an absolute position
            # more than 400 mm from the origin
            if !ismissing(R_L) && abs(R_L) > 100
                flag += sign(R_L)
            elseif (Rval > 400) === true
                flag += 1
            elseif (Lval > 400) === true
                flag -= 1
            end
        end

        # Code indeterminate cases as '0', to be reconsidered later
        side[i] = iszero(flag) ? '0' :
                     (flag > 0 ? 'R' : 'L')
    end

    # Detect and remove double/multiple taps on the touch target (debounce, ish)
    cycle_times = diff(pks)
    med_length = median(cycle_times)
    if any(<(med_length/2), cycle_times)
        shorti = findall(<(med_length/2), cycle_times)

        # The second touch in a double touch will be the "bad" touch, so iterate in reverse
        # order to avoid changing future iterates
        for i in reverse(shorti)
            if side[i+1] == side[i]
                @debug "Double touch at $(pks[i+1])" trial
                deleteat!(pks, i+1)
                deleteat!(side, i+1)
            end
        end
    end

    cycle_times = diff(pks)
    med_length = median(cycle_times)

    # Investigate long cycles (which could be due to missed touches, and reaches could
    # recovered using kinematic data). This stage is identical (and therefore unhelpful) for
    # the trial with bad touch signal where the existing events (at this point) were created
    # with the same finger marker data
    longi = findall(>(med_length*1.5), cycle_times)
    if !isempty(longi)
        rmkr = @view rsrc.point["RFIN"][:,2]
        lmkr = @view rsrc.point["LFIN"][:,2]
        ridxs = argmaxima(rmkr, 5; strict=false);
        peakproms!(ridxs, rmkr; strict=false, minprom=100)
        lidxs = argmaxima(lmkr, 5; strict=false);
        peakproms!(lidxs, lmkr; strict=false, minprom=100)
        kin_pks = sort!(convert(Vector{Float32}, [ridxs; lidxs]))

        # Match new kinematically determined events with existing events (which will be from
        # the touch target in most trials). Since there will be differences in the timing of
        # touch events vs kinematic events, record the average error offset between the two
        # types of events and remove from any newly added kinematic events (to avoid
        # introducing timing inconsistencies from determining events with different methods)
        _, _, errors, _ = matchevents(kin_pks, convert(Vector{Float32}, pks))
        mean_error = mean(errors)
        for i in reverse(longi) # Iterate in reverse order again to avoid changing future iterates
            new_pks = filter(x -> pks[i]+med_length/2 < x < pks[i+1]-med_length/2, kin_pks)
            if !isempty(new_pks)
                new_pks .+= mean_error
                # @debug "Adding pks at $(round.(new_pks./10; digits=1)))"
                splice!(pks, i+1:0, round.(Int, new_pks))
            end
        end
    end

    cycle_times = diff(pks)
    med_length = median(cycle_times)
    # Stopping reaching more than 5 sec before the end of the trial is suspicious
    if last(pks) < numpointframes(rsrc) - 500
        # Use all right and left markers that would show clear peaks from reaching, and then
        # cluster peak timing (with `dbscan`) for greater confidence in individual timing of
        # peaks and to reduce differences in peak timing between markers
        present_mkrs = find_present_markers(rsrc,
            ["LFIN", "LIDX", "LWRM", "LWRL", "RFIN", "RIDX", "RWRM", "RWRL"];
            threshold=0.25)
        mkr_peaks = Dict()
        foreach(present_mkrs) do mkrname
            _pks = get!(mkr_peaks, mkrname) do
                @views argmaxima(rsrc.point[mkrname][:,2], 5; strict=false)
            end
            # @show "before" length(_pks)
            @views peakproms!(_pks, rsrc.point[mkrname][:,2]; strict=false, minprom=100)
            # @show "after" length(_pks)
        end

        # We only care about "new" peaks, ie peaks that occur during long cycles where we
        # suspect there may be undiscovered reaches
        lpks = sort!(reduce(vcat, (mkr_peaks[mkrname]
            for mkrname in filter(contains(r"^L"), present_mkrs))))
        filter!(>((last(pks)+med_length/2)), lpks)
        rpks = sort!(reduce(vcat, (mkr_peaks[mkrname]
            for mkrname in filter(contains(r"^R"), present_mkrs))))
        filter!(>((last(pks)+med_length/2)), rpks)

        # Cluster any peaks that occur within 100 ms
        lclust = dbscan(reshape(Float32.(lpks), 1, :), 10; min_cluster_size=2)
        rclust = dbscan(reshape(Float32.(rpks), 1, :), 10; min_cluster_size=2)

        new_pks = Int[]
        if !isempty(lclust.clusters)
            # Round new peaks from clusters to the nearest frame and add to array of
            # confirmed new peaks
            act_lpks = map(lclust.clusters) do clust
                round(Int, mean(lpks[clust.core_indices]))
            end
            # Remove peaks in (successful) clusters from the array of candidate "new" peaks
            # (i.e. `lpks`)
            deleteat!(lpks, mapreduce(x -> x.core_indices, vcat, lclust.clusters; init=Int[]))

            # Remove any peaks in the array of "candidate" new peaks (i.e. `lpks`) that are
            # too close to any clusters (i.e. peaks that were outside the cluster radius of
            # 100 ms, but otherwise too near to plausibly be a correct/actual reach)
            ltree = KDTree(reshape(Float32.(lpks), 1, :))
            tooclose = sort!(unique!(reduce(vcat, inrange(ltree, reshape(Float32.(act_lpks), 1, :), 15))))
            !isempty(tooclose) && deleteat!(lpks, tooclose)

            # Append new peaks from the clusters
            append!(new_pks, act_lpks)
        end
        if !isempty(rclust.clusters) # Same as above, but for the right side
            act_rpks = map(rclust.clusters) do clust
                round(Int, mean(rpks[clust.core_indices]))
            end

            deleteat!(rpks, mapreduce(x -> x.core_indices, vcat, rclust.clusters; init=Int[]))
            rtree = KDTree(reshape(Float32.(rpks), 1, :))
            tooclose = sort!(unique!(reduce(vcat, inrange(rtree, reshape(Float32.(act_rpks), 1, :), 15))))
            !isempty(tooclose) && deleteat!(rpks, tooclose)
            append!(new_pks, act_rpks)
        end

        # Use any remaining single kinematic peaks (i.e. peak from a single marker that
        # did not match to a cluster, and occured with reasonable temporal separation from
        # any other extant events; anything that survived the above pruning).
        append!(new_pks, lpks, rpks)

        sort!(new_pks)
        # @assert minimum(diff(pks)) > 15
        # println("unclust = ", new_pks)
        append!(pks, new_pks) # These new peaks will only be part of the last 5 sec of the trial
    end

    cycle_times = diff(pks)
    med_length = median(cycle_times)

    # Re-evaluate side (same alg as previous) since more peaks have been added
    resize!(side, length(pks))
    for (i, touch) in enumerate(pks)
        flag = 0
        for (R, L) in zip(["RIDX", "RFIN", "RWRM", "RWRL"], ["LIDX", "LFIN", "LWRM", "LWRL"])
            Rval = get(rsrc.point, R, LazyFill(missing))[max(touch, 1), 2]
            Lval = get(rsrc.point, L, LazyFill(missing))[max(touch, 1), 2]
            R_L = Rval - Lval
            if !ismissing(R_L)# && abs(R_L) > 0
                flag += sign(R_L)
            elseif (Rval > 400) === true
                flag += 1
            elseif (Lval > 400) === true
                flag -= 1
            end
        end
        side[i] = iszero(flag) ? '0' :
                     (flag > 0 ? 'R' : 'L')
    end

    # Detect and remove double/multiple pks
    intvls = intervals(pks)
    med_length = median(length.(intvls))

    shorti = findall(<(med_length)∘length, intvls)
    if !isempty(shorti)
        present_mkrs = find_present_markers(rsrc,
            ["LFIN", "LIDX", "LWRM", "LWRL", "RFIN", "RIDX", "RWRM", "RWRL"];
            threshold=0.25)
        # Different threshold this time because any new short cycles are likely due to the
        # kinematic touches, and may not be actual double touches (i.e. the 2 data methods
        # disagree enough to appear as separate and repeated touches)
        for i in reverse(shorti) # Delete from the back to avoid changing earlier indices
            if side[i+1] == side[i] # Double touch
                # For each hand/reach event, calculate the average AP position from all
                # markers present (i.e. non-missing) at that moment
                inline_switch(f, val, init) = f(val) ? val : init
                l_mkrsi = mean(inline_switch(!isempty, skipmissing([ rsrc.point[mkrname][pks[i], 2]
                    for mkrname in filter(contains(r"^L"), present_mkrs) ]), [missing]))
                r_mkrsi = mean(inline_switch(!isempty, skipmissing([ rsrc.point[mkrname][pks[i], 2]
                    for mkrname in filter(contains(r"^R"), present_mkrs) ]), [missing]))
                l_mkrsi1 = mean(inline_switch(!isempty, skipmissing([ rsrc.point[mkrname][pks[i+1], 2]
                    for mkrname in filter(contains(r"^L"), present_mkrs) ]), [missing]))
                r_mkrsi1 = mean(inline_switch(!isempty, skipmissing([ rsrc.point[mkrname][pks[i+1], 2]
                    for mkrname in filter(contains(r"^R"), present_mkrs) ]), [missing]))

                # Noisy markers and/or gaps could create peaks with prominence > 100 that
                # aren't reflective of an actual reach
                if side[i] == 'R'
                    reachi = r_mkrsi - l_mkrsi
                    reachi1 = r_mkrsi1 - l_mkrsi1
                    # If there's more than a 100 mm AP distance between hand positions
                    # between the two reaches, the shorter reach should be removed
                    if coalesce(abs(reachi - reachi1), 0) > 100
                        # Remove smaller reach
                        j = argmin((reachi, reachi1))-1
                        @debug "Double touch (false peak) at $(pks[i+j])" trial
                        deleteat!(pks, i+j)
                        deleteat!(side, i+j)
                    end
                elseif side[i] == 'L'
                    reachi = l_mkrsi - r_mkrsi
                    reachi1 = l_mkrsi1 - r_mkrsi1
                    if coalesce(abs(reachi - reachi1), 0) > 100
                        j = argmin((reachi, reachi1))-1
                        @debug "Double touch (false peak) at $(pks[i+j])" trial
                        deleteat!(pks, i+j)
                        deleteat!(side, i+j)
                    end
                else
                    # Some (few) subjects occasionally did complete/actual double touches
                    # (e.g. due to missing the touch target, and re-reaching to ensure a
                    # percieved "correct" reach). Remove the second touch
                    @debug "Double touch (kin) at $(pks[i+1])" trial
                    deleteat!(pks, i+1)
                    deleteat!(side, i+1)
                end
            end
        end
    end

    # All alternatives to discover new reaches have been used, therefore, we assume that all
    # reaches have in fact been discovered, and any reaches where the side is ambiguous are
    # filled with whichever side fits for alternating reaches
    if any(==('0'), side)
        z0s = findall(==('0'), side)
        for i in z0s
            i == 1 && continue
            side[i] = (side[i-1] == 'R') ? 'L' : 'R' # Opposite side from last reach
        end

        rz0s = findall(==('0'), side)
        for i in rz0s
            i == lastindex(side) && continue
            side[i] = (side[i+1] == 'R') ? 'L' : 'R'
        end
    end

    for i in reverse(axes(side, 1)[2:end])
        if side[i] == side[i-1]
            # Double touches at this point are actual double touches
            @debug "Double touch (final) on $(side[i]) at $(pks[i])" trial
            deleteat!(pks, i)
            deleteat!(side, i)
        end
    end

    # Assert no remaining unknown reaches
    @assert !any(==('0'), side)

    # Reaching is supposed to begin at 5 sec into trial; ignore/remove earlier reaches
    if filter5
        early = findall(<(500), pks)
        deleteat!(pks, early)
        deleteat!(side, early)
    end

    intvls = intervals(pks)
    med_length = median(length.(intvls))

    # Final reach intervals must be complete "stride" equivalents (e.g. two reaches for a
    # given side, separated by a reach by the other side), where the intervals between both
    # "step" equivalents are less than the "plausible freezing episode" threshold
    left_reaches = CycleInterval{UnitRange{Int}}[]
    right_reaches = CycleInterval{UnitRange{Int}}[]
    for i in eachindex(pks, side)[2:end-1]
        if side[i-1] == 'L' && side[i] == 'R' && side[i+1] == 'L'
            if pks[i] - pks[i-1] ≤ med_length*2 && pks[i+1] - pks[i] ≤ med_length*2
                push!(right_reaches, CycleInterval(pks[i-1], pks[i], pks[i+1]))
            end
        elseif side[i-1] == 'R' && side[i] == 'L' && side[i+1] == 'R'
            if pks[i] - pks[i-1] ≤ med_length*2 && pks[i+1] - pks[i] ≤ med_length*2
                push!(left_reaches, CycleInterval(pks[i-1], pks[i], pks[i+1]))
            end
        end
    end

    return left_reaches, right_reaches, any(>(med_length*2), length.(intvls))
end


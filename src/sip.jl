function condition_FP(sig)
    df = let
        n = 4
        bw = Butterworth(n)
        # Cutoff at 50Hz to maintain sharp corner of foot-strike and liftoff
        lpf = Lowpass(50.; fs=1000)
        digitalfilter(lpf, bw)
    end

    return filtfilt(df, sig)
end

"Calculate velocity magnitude (i.e. norm'ed) of a marker"
function marker_velocity(marker; dt=1)
    diffed = centraldiff(marker; dt)
    vel = Vector{Union{Missing,Float32}}(missing, size(diffed, 1))
    for i in axes(diffed, 1)
        if ismissing(diffed[i,1])
            vel[i] = missing
        else
            buf = SVector{3,Float32}(diffed[i,1]::Float32, diffed[i,2]::Float32,
                diffed[i,3]::Float32)
            vel[i] = norm(buf)
        end
    end
    return vel
end

"""
    markergroup_velocity(markers...; dt=1)

Calculate the average absolute velocity of a group of markers, progressively allowing for
any/multiple markers to be missing at any instant.

(i.e. Average velocity magnitude of all markers in group that are present at any given
instant)
"""
function markergroup_velocity(markers...; dt=1)
    diffs = centraldiff.(markers; dt)

    mdiffs = fill(NaN32, size(diffs[1]))
    counts = zeros(Int, size(diffs[1], 1))
    for i in eachindex(diffs)
        nmissi = findall(!ismissing, view(diffs[i], :, 1))
        mdiffs[nmissi,:] .= 0 # remove poison
        counts[nmissi] .+= 1 # count for future averaging
    end
    for i in eachindex(diffs)
        nmissi = findall(!ismissing, view(diffs[i], :, 1))
        mdiffs[nmissi,:] .+= diffs[i][nmissi,:] # sum all present
    end
    mdiffs ./= counts

    vel = Vector{Union{Missing,Float32}}(missing, size(mdiffs, 1))
    for i in axes(mdiffs, 1)
        if ismissing(mdiffs[i,1])
            vel[i] = missing
        else
            buf = SVector{3,Float32}(mdiffs[i,1]::Float32, mdiffs[i,2]::Float32,
                mdiffs[i,3]::Float32)
            vel[i] = norm(buf)
        end
    end
    return vel
end

"""
    findmarker_stanceheight(mkr)

Calculate the height of a (typically foot) marker when the foot is fully planted and
not-moving, identified with a functional 10%ile velocity threshold.
"""
function findmarker_stanceheight(mkr)
    # Function is only called on markers that are present at least 50% of the trial
    # The markers present 50-55% are iffy enough to skip the functional threshold and just
    # use the minimum vertical marker position.
    if 0.5 ≤ count(!ismissing, @view(mkr[:,1]))/size(mkr, 1) ≤ 0.55
        return minimum(skipmissing(@view(mkr[:,3])))
    else
        vel = running(mean, marker_velocity(mkr; dt=inv(100)), 10)
        # Find all instances when marker is present with a velocity less than the 10%ile
        q10 = quantile(skipmissing(vel), 0.1)
        is = findall(x -> (x < q10) === true, vel)
        return median(@view(mkr[is,3]))
    end
end

"""
    find_present_markers(c, mkrs; threshold=0.5)

Find markers from `mkrs` in the C3D file `c` that are present at least `threshold`% of the
trial
"""
function find_present_markers(c, mkrs; threshold=0.5)
    present_mkrs = filter(x -> haskey(c.point, x), mkrs)
    filter!(present_mkrs) do mkr
        count(!ismissing, @view(c.point[mkr][:,1]))/size(c.point[mkr], 1) > threshold
    end

    present_mkrs
end

# Not used in find_cycles(trial, src, Val(:sip))
"Remove the stance height from the present markers"
function zero_present_markers(c, mkrs)
    present_mkrs = find_present_markers(c, mkrs)
    return zero_markers(c, present_mkrs)
end

"Remove the stance height from `mkrs`"
function zero_markers(c, mkrs)
    zeroed_mkrs = Dict()
    foreach(mkrs) do mkr
        get!(zeroed_mkrs, mkr) do
            c.point[mkr][:,3] .- findmarker_stanceheight(c.point[mkr])
        end
    end

    return zeroed_mkrs
end

function find_cycles(trial, src, ::Val{:sip}; filter5=true, peak_separation=15, min_stepheight=8)
    rsrc = readsource(src; strip_prefixes=true)

    # Determine initial/base events from feet and ankle markers
    present_footmkrs = filter(x -> haskey(rsrc.point, x),
        ["LANK", "LHEE", "LMT5", "LTOE", "RANK", "RHEE", "RMT5", "RTOE"])
    filter!(present_footmkrs) do mkr
        count(!ismissing, @view(rsrc.point[mkr][:,1]))/size(rsrc.point[mkr], 1) > 0.5
    end

    mkr_height = Dict()
    foreach(present_footmkrs) do mkrname
        get!(mkr_height, mkrname) do
            findmarker_stanceheight(rsrc.point[mkrname])
        end
    end

    mkr_peaks = Dict()
    foreach(present_footmkrs) do mkrname
        pks = get!(mkr_peaks, mkrname) do
            @views argmaxima(rsrc.point[mkrname][:,3], peak_separation; strict=false)
        end
        peakheights!(pks, rsrc.point[mkrname][pks, 3];
            minheight=mkr_height[mkrname]+min_stepheight)
    end

    lpks = sort!(reduce(vcat, (mkr_peaks[mkrname]
         for mkrname in filter(contains(r"^L"), present_footmkrs))));
    rpks = sort!(reduce(vcat, (mkr_peaks[mkrname]
        for mkrname in filter(contains(r"^R"), present_footmkrs))));

    # For peaks from the same side, form clusters of at least 2 events that occur within
    # 150 ms of each other
    lclust = dbscan(reshape(Float32.(lpks), 1, :), peak_separation; min_cluster_size=2)
    rclust = dbscan(reshape(Float32.(rpks), 1, :), peak_separation; min_cluster_size=2)

    # Round center (time) of clusters to the nearest frame (aka Integer)
    act_lpks = map(lclust.clusters) do clust
        round(Int, mean(lpks[clust.core_indices]))
    end
    act_rpks = map(rclust.clusters) do clust
        round(Int, mean(rpks[clust.core_indices]))
    end

    # Create array of confirmed steps and an array for the side of each step
    pks = [act_lpks; act_rpks]
    ix = sortperm(pks)
    permute!(pks, ix)
    side = permute!([ fill('L', length(act_lpks)); fill('R', length(act_rpks)) ], ix)

    # Investigate long cycles (which could be due to missed steps)
    intvls = intervals(pks)
    med_length = median(length.(intvls))
    longi = findall(>(med_length*1.5)∘length, intvls)
    if !isempty(longi)
        # Remove peaks from initial arrays that were clustered and added to `pks`
        filt_lpks = deleteat!(copy(lpks), mapreduce(x -> x.core_indices, vcat, lclust.clusters))
        filt_rpks = deleteat!(copy(rpks), mapreduce(x -> x.core_indices, vcat, rclust.clusters))

        # SUMMARY:
        # For all long cycles:
        # 1. Check for remaining acceptably high peaks (ie above the `min_stepheight` in
        # lpks/rpks) that would fill the long cycles
            # - Must be more than 1/4 of the median cycle length away from the next/prev
            # event to add
            # - No "double steps" (ie within 150 ms) should exist since any group of >2
            # would be clustered already
        # 2. Check for peaks that are within the
            # respective FP's and show appropriate loading/unloading (defined as difference
            # of <5%BW and >95%BW for the stepping and stance feet, respectively)

        # First pass check
        for i in reverse(longi)
            # cddt = candidate
            lcddt = intersect(first(intvls[i])+med_length/4:last(intvls[i])-med_length/4, filt_lpks)
            rcddt = intersect(first(intvls[i])+med_length/4:last(intvls[i])-med_length/4, filt_rpks)
            setdiff!(filt_lpks, lcddt)
            setdiff!(filt_rpks, rcddt)

            cddts = [lcddt; rcddt]
            ix = sortperm(cddts)
            permute!(cddts, ix)
            new_sides = permute!([ fill('L', length(lcddt)); fill('R', length(rcddt)) ], ix)

            splice!(pks, (i+1):i, cddts)
            splice!(side, (i+1):i, new_sides)
        end
        @assert issorted(pks)

        # Second pass checks
        RFP_X_mn = rsrc.groups[:FORCE_PLATFORM][:CORNERS][1,1,1]
        RFP_X_mx = rsrc.groups[:FORCE_PLATFORM][:CORNERS][1,2,1]
        RFP_Y_mn = rsrc.groups[:FORCE_PLATFORM][:CORNERS][2,2,1]
        RFP_Y_mx = rsrc.groups[:FORCE_PLATFORM][:CORNERS][2,3,1]
        LFP_X_mn = rsrc.groups[:FORCE_PLATFORM][:CORNERS][1,1,2]
        LFP_X_mx = rsrc.groups[:FORCE_PLATFORM][:CORNERS][1,2,2]
        LFP_Y_mn = rsrc.groups[:FORCE_PLATFORM][:CORNERS][2,2,2]
        LFP_Y_mx = rsrc.groups[:FORCE_PLATFORM][:CORNERS][2,3,2]

        # Reduce size of FPs by 5.5 cm for shoe clearance (assumption: markers are radially
        # a max of 5.5 cm from shoe edge)
        RFP_X_mn += copysign(50, RFP_X_mx - RFP_X_mn)
        RFP_X_mx -= copysign(50, RFP_X_mx - RFP_X_mn)
        RFP_Y_mn += copysign(50, RFP_Y_mx - RFP_Y_mn)
        RFP_Y_mx -= copysign(50, RFP_Y_mx - RFP_Y_mn)
        LFP_X_mn += copysign(50, LFP_X_mx - LFP_X_mn)
        LFP_X_mx -= copysign(50, LFP_X_mx - LFP_X_mn)
        LFP_Y_mn += copysign(50, LFP_Y_mx - LFP_Y_mn)
        LFP_Y_mx -= copysign(50, LFP_Y_mx - LFP_Y_mn)

        RFP_X_mn, RFP_X_mx = minmax(RFP_X_mn, RFP_X_mx)
        RFP_Y_mn, RFP_Y_mx = minmax(RFP_Y_mn, RFP_Y_mx)
        LFP_X_mn, LFP_X_mx = minmax(LFP_X_mn, LFP_X_mx)
        LFP_Y_mn, LFP_Y_mx = minmax(LFP_Y_mn, LFP_Y_mx)

        RFP = rsrc.analog["Force.Fz1"]
        LFP = rsrc.analog["Force.Fz2"]

        usedL = similar(filt_lpks, 0)
        for i in eachindex(filt_lpks)
            if -RFP[filt_lpks[i]*10] ≥ conditions(trial)[:bodyweight]*9.81*.9 && # 90%BW Load
                -LFP[filt_lpks[i]*10] ≤ conditions(trial)[:bodyweight]*9.81*.05 # 5%BW Load
                rmkrs = [ rsrc.point[mkr][filt_lpks[i],:] for mkr in filter(contains(r"^R"), present_footmkrs) ]

                if all(x -> RFP_X_mn ≤ coalesce(x[1], RFP_X_mn) ≤ RFP_X_mx && RFP_Y_mn ≤ coalesce(x[2], RFP_Y_mn) ≤ RFP_Y_mx, rmkrs)
                    j = searchsortedfirst(pks, filt_lpks[i])
                    insert!(pks, j, filt_lpks[i])
                    insert!(side, j, 'L')
                    push!(usedL, i)
                else
                end
            end
        end
        !isempty(usedL) && deleteat!(filt_lpks, usedL)
        usedR = similar(filt_lpks, 0)
        for i in eachindex(filt_rpks)
            if -LFP[filt_rpks[i]*10] ≥ conditions(trial)[:bodyweight]*9.81*.9 && # 90%BW Load
                -RFP[filt_rpks[i]*10] ≤ conditions(trial)[:bodyweight]*9.81*.05 # 5%BW Load
                lmkrs = [ rsrc.point[mkr][filt_rpks[i],:] for mkr in filter(contains(r"^L"), present_footmkrs) ]

                if all(x -> LFP_X_mn ≤ coalesce(x[1], LFP_X_mn) ≤ LFP_X_mx && LFP_Y_mn ≤ coalesce(x[2], LFP_Y_mn) ≤ LFP_Y_mx, lmkrs)
                    j = searchsortedfirst(pks, filt_rpks[i])
                    insert!(pks, j, filt_rpks[i])
                    insert!(side, j, 'L')
                    push!(usedR, i)
                else
                end
            end
        end
        !isempty(usedR) && deleteat!(filt_rpks, usedR)
    end

    # Detect and remove any potential double steps
    intvls = intervals(pks)
    med_length = median(length.(intvls))
    shorti = findall(<(med_length)∘length, intvls)
    if !isempty(shorti)
        for i in reverse(shorti)
            if side[i+1] == side[i] # Double step detected, context not established
                #    ⬇ check that i-1 is a valid index since i is decreasing (reverse iteration)
                if i > 1 && length(intvls[i-1]) < med_length # Double short step
                    @debug "REMOVED: (Early) double step at $(pks[i])" trial
                    deleteat!(pks, i)
                    deleteat!(side, i)
                else
                    # Double step, normal length; test for incongruent marker positions
                    _mkrheights = [ rsrc.point[mkrname][pks[i+1], 3] - mkr_height[mkrname]
                        for mkrname in present_footmkrs ]
                    if count(<(5), skipmissing(_mkrheights)) > 2
                        # At least 2 markers are less than 5mm above their stance height;
                        # foot is not actually raised
                        @debug "REMOVED: Double step at $(pks[i+1])" trial
                        deleteat!(pks, i+1)
                        deleteat!(side, i+1)
                    else
                        # Indeterminate double step; warn
                        @debug "Double step at $(pks[i+1])" trial
                    end
                end
            elseif length(intvls[i]) ≤ 20 # Not a double (same-side) step, but too short to be good
                # Calculate zeroed heights (ie subtract the stance height)
                l_mkrheightsi = mean(skipmissing([ rsrc.point[mkrname][pks[i], 3] - mkr_height[mkrname]
                    for mkrname in filter(contains(r"^L"), present_footmkrs) ]))
                r_mkrheightsi = mean(skipmissing([ rsrc.point[mkrname][pks[i], 3] - mkr_height[mkrname]
                    for mkrname in filter(contains(r"^R"), present_footmkrs) ]))
                l_mkrheightsi1 = mean(skipmissing([ rsrc.point[mkrname][pks[i+1], 3] - mkr_height[mkrname]
                    for mkrname in filter(contains(r"^L"), present_footmkrs) ]))
                r_mkrheightsi1 = mean(skipmissing([ rsrc.point[mkrname][pks[i+1], 3] - mkr_height[mkrname]
                    for mkrname in filter(contains(r"^R"), present_footmkrs) ]))

                # If one side is higher than the other side for both "steps" than the lower
                # "step" is a false-positive and should be removed
                if l_mkrheightsi > r_mkrheightsi && l_mkrheightsi1 > r_mkrheightsi1
                    j = findnext(==('R'), side, i)
                    @assert j == i || j == i+1
                    deleteat!(pks, j)
                    deleteat!(side, j)
                elseif r_mkrheightsi > l_mkrheightsi && r_mkrheightsi1 > l_mkrheightsi1
                    j = findnext(==('L'), side, i)
                    @assert j == i || j == i+1
                    deleteat!(pks, j)
                    deleteat!(side, j)
                end
            end
        end
    end

    # Stepping is supposed to begin at 5 sec into trial; ignore/remove earlier steps
    if filter5
        early = findall(<(500), pks)
        deleteat!(pks, early)
        deleteat!(side, early)
    end

    intvls = intervals(pks)
    med_length = median(length.(intvls))

    # Final steps intervals must be complete strides, where both step times are less than
    # the "plausible freezing episode" threshold
    left_steps = CycleInterval{UnitRange{Int}}[]
    right_steps = CycleInterval{UnitRange{Int}}[]
    for i in eachindex(pks, side)[2:end-1]
        if side[i-1] == 'L' && side[i] == 'R' && side[i+1] == 'L'
            if pks[i] - pks[i-1] ≤ med_length*1.75 && pks[i+1] - pks[i] ≤ med_length*1.75
                push!(right_steps, CycleInterval(pks[i-1], pks[i], pks[i+1]))
            end
        elseif side[i-1] == 'R' && side[i] == 'L' && side[i+1] == 'R'
            if pks[i] - pks[i-1] ≤ med_length*1.75 && pks[i+1] - pks[i] ≤ med_length*1.75
                push!(left_steps, CycleInterval(pks[i-1], pks[i], pks[i+1]))
            end
        end
    end

    return left_steps, right_steps, any(>(med_length*1.75), length.(intvls))
end


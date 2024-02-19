module PaperResearchCode

using DatasetManager, C3D, DSP, Peaks, Statistics, StatsBase, Biomechanics, StaticArrays,
    LinearAlgebra, RollingFunctions, NearestNeighbors, Clustering

export condition_FP, condition_emg, condition_marker, condition_touch, find_cycles,
    marker_velocity, findmarker_stanceheight, CycleInterval, midpoint, extents,
    contiguousranges, invert_array_of_ranges, zero_present_markers, find_present_markers,
    zero_markers, cycles_analysis

include("general.jl")
include("rrt.jl")
include("sip.jl")
if isfile("emg.jl")
    include("emg.jl")
end

end

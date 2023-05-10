"""
    track(
    imgsdir::String,
    propsdir::String,
    timedeltasdir::String,
    paramsdir::String,
    outdir::String,
)

$(include("track-docstring.jl"))
"""
function track(;
    imgs::String, props::String, deltat::String, params::String, output::String, args...
)
    @show imgs props deltat params output
    # imgs = deserialize(joinpath(imgsdir, "segmented_floes.jls"))
    # properties = deserialize(joinpath(propsdir, "floe_props.jls"))
    # delta_time = deserialize(joinpath(timedeltasdir, "passtimes.jls")) # TODO: process passtimes
    # params = TOML.parsefile(joinpath(paramsdir, "tracker-params.toml"))
    # t1, t2, t3, mc_thresholds =
    #     dict2nt.(params["t1"], params["t2"], params["t3"], params["mc_thresholds"])
    # condition_thresholds = (t1=t1, t2=t2, t3=t3)
    # pairs = pairfloes(imgs, properties, delta_time, condition_thresholds, mc_thresholds)
    # serialize(joinpath(outdir, "matched_pairs.jls"), pairs)
    return nothing
end

"""
    dict2nt(d)

Convert a dictionary `d` to a NamedTuple.
"""
dict2nt(d) = NamedTuple((Symbol(key), value) for (key, value) in d)

function parse_item(::Type{Vector{Int64}}, x::AbstractString)
    return parse.(Int64, split(x))
end

function parse_item(::Type{Vector{Float64}}, x::AbstractString)
    return parse.(Float64, split(x))
end
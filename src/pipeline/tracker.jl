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
function track(
    imgsdir::String,
    propsdir::String,
    timedeltasdir::String,
    paramsdir::String,
    outdir::String,
)
    imgs = deserialize(joinpath(imgsdir, "segmented_floes.jls"))
    properties = deserialize(joinpath(propsdir, "floe_props.jls"))
    delta_time = deserialize(joinpath(timedeltasdir, "passtimes.jls")) # TODO: process passtimes
    params = TOML.parsefile(joinpath(paramsdir, "tracker-params.toml"))
    t1, t2, t3, mc_thresholds =
        dict2nt.(params["t1"], params["t2"], params["t3"], params["mc_thresholds"])
    condition_thresholds = (t1=t1, t2=t2, t3=t3)
    pairs = pairfloes(imgs, properties, delta_time, condition_thresholds, mc_thresholds)
    serialize(joinpath(outdir, "matched_pairs.jls"), pairs)
    return nothing
end

"""
    dict2nt(d)

Convert a dictionary `d` to a NamedTuple.
"""
dict2nt(d) = NamedTuple((Symbol(key), value) for (key, value) in d)

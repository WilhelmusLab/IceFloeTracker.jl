"""
    pair_floes(
    imgspath::String, propspath::String, condition_thresholds, mc_thresholds=(area3=0.18, area2=0.236, corr=0.68)
)

$(include("pair_floes_docstring.jl"))

"""
function pair_floes(
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
    t1, t2, t3 = dict2nt.(params["t1"], params["t2"], params["t3"])
    condition_thresholds = (t1=t1, t2=t2, t3=t3)
    mc_thresholds = dict2nt(params["mc_thresholds"])
    pairs = pairfloes(imgs, properties, delta_time, condition_thresholds, mc_thresholds)
    serialize(joinpath(outdir, "matched_pairs.jls"), pairs)
    return nothing
end

"""
    dict2nt(d)

Converts a dictionary `d` to a NamedTuple.
"""
dict2nt(d) = NamedTuple((Symbol(key), value) for (key, value) in d)
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
function track(; args...)
    condition_thresholds, mc_thresholds = get_thresholds(; args...)
    vals = values(args)
    imgs = deserialize(joinpath(vals.imgs, "segmented_floes.jls"))
    props = deserialize(joinpath(vals.props, "floe_props.jls"))
    deltat = deserialize(joinpath(vals.deltat, "timedeltas.jls"))
    @show vals.output
    serialize(
        joinpath(vals.output, "tracked_floes.jls"),
        IceFloeTracker.pairfloes(imgs, props, deltat, condition_thresholds, mc_thresholds),
    )
    return nothing
end

function parse_params(params::AbstractString)
    params = parsefile(params)
    area = params["area"]
    t1 = dict2nt(params["t1"])
    t2 = (area=area, dict2nt(params["t2"])...)
    t3 = (area=area, dict2nt(params["t3"])...)
    condition_thresholds = (t1=t1, t2=t2, t3=t3)
    d = dict2nt(params["mc_thresholds"])
    @show d
    mc_thresholds = mkmct(d)
    return condition_thresholds, mc_thresholds
end

function parse_params(; args...)
    d = values(args)
    condition_thresholds = (
        t1=(dt=parselistas(Int64, d.dt_thresh), dist=parselistas(Int64, d.dist)),
        t2=(
            area=d.area,
            arearatio=d.Larearatio,
            convexarearatio=d.Lconvexarearatio,
            majaxisratio=d.Lmajaxisratio,
            minaxisratio=d.Lminaxisratio,
        ),
        t3=(
            area=d.area,
            arearatio=d.Sarearatio,
            convexarearatio=d.Sconvexarearatio,
            majaxisratio=d.Smajaxisratio,
            minaxisratio=d.Sminaxisratio,
        ),
    )
    mc_thresholds = mkmct(d)
    return condition_thresholds, mc_thresholds
end

function mkmct(d)
    return (
        comp=(mxrot=d.mxrot, sz=d.sz, comp=d.comp, mm=d.mm, psi=d.psi),
        goodness=(corr=d.corr, area2=d.area2, area3=d.area3),
    )
end

function get_thresholds(; args...)
    v = values(args)
    if !isnothing(v.params)
        condition_thresholds, mc_thresholds = parse_params(v.params)
    else
        condition_thresholds, mc_thresholds = parse_params(; args...)
    end
    return condition_thresholds, mc_thresholds
end

"""
    dict2nt(d)

Convert a dictionary `d` to a NamedTuple.
"""
dict2nt(d) = NamedTuple((Symbol(key), value) for (key, value) in d)

parselistas(T, x) = [parse(T, i) for i in split(x)]

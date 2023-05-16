include("../scripts/track-cli-config.jl")
settings = ArgParseSettings(; autofix_names=true)
mkclitrack!(settings)

params_path = "test_inputs/pipeline/tracker/tracker-params.toml"
data_path = "test_inputs/tracker/tracker_test_data.dat"
temp = mkpath(joinpath(@__DIR__, "__temp__"))
cliparams = [
    "track",
    "--imgs",
    temp,
    "--props",
    temp,
    "--deltat",
    temp,
    "--output",
    temp,
    "--params",
    params_path,
]

obj = deserialize(data_path)

serialize(joinpath(temp, "segmented_floes.jls"), obj.imgs)
serialize(joinpath(temp, "floe_props.jls"), obj.props)
serialize(joinpath(temp, "timedeltas.jls"), [15, 20])

argsparam = parse_args(cliparams, settings; as_symbols=true)
cmd = argsparam[:_COMMAND_]
argsparam[cmd]
IceFloeTracker.Pipeline.track(; argsparam[cmd]...)
tracked = deserialize(joinpath(temp, "tracked_floes.jls"))
@test isfile(joinpath(temp, "tracked_floes.jls"))
@test length(tracked) == 2
@test nrow(tracked[1].props1) == 9
@test nrow(tracked[1].props2) == 9
@test nrow(tracked[2].props1) == 21
@test nrow(tracked[2].props2) == 21
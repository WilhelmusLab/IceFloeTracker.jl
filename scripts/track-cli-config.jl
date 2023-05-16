"""
    mkclitrack!(settings)

Set up command line interface for the `track` command.
"""
function mkclitrack!(settings)
    @add_arg_table! settings begin
        "track"
        help = "Pair ice floes in day k with ice floes in day k+1"
        action = :command
    end

    add_arg_group!(settings["track"], "arguments")
    @add_arg_table! settings["track"] begin
        "--imgs"
        help = "Path to object with segmented images"
        required = true

        "--props"
        help = "Path to object with extracted features"
        required = true

        "--deltat"
        help = "Path to object with time deltas between consecutive images"
        required = true

        "--output", "-o"
        help = "Output directory"
        required = true
    end

    add_arg_group!(settings["track"], "optional arguments")
    @add_arg_table! settings["track"] begin
        "--params"
        help = "Path to TOML file with algorithm parameters"

        "--area"
        help = "Area thresholds to use for pairing floes"
        arg_type = Int64
        default = 1200

        "--dist"
        help = "Distance threholds to use for pairing floes"
        default = "15 30 120"

        "--dt-thresh"
        help = "Time thresholds to use for pairing floes"
        default = "30 100 1300"

        "--Sarearatio"
        help = "Area ratio threshold"
        arg_type = Float64
        default = 0.18

        "--Smajaxisratio"
        help = "Major axis ratio threshold to use for pairing floes"
        arg_type = Float64
        default = 0.1

        "--Sminaxisratio"
        help = "Minor axis ratio threshold to use for pairing floes"
        arg_type = Float64
        default = 0.12

        "--Sconvexarearatio"
        help = "Convex area ratio threshold to use for pairing floes"
        arg_type = Float64
        default = 0.14

        "Larearatio"
        help = "Area ratio threshold"
        arg_type = Float64
        default = 0.28

        "--Lmajaxisratio"
        help = "Major axis ratio threshold to use for pairing floes"
        arg_type = Float64
        default = 0.1

        "--Lminaxisratio"
        help = "Minor axis ratio threshold to use for pairing floes"
        arg_type = Float64
        default = 0.15

        "--Lconvexarearatio"
        help = "Convex area ratio threshold to use for pairing floes"
        arg_type = Float64
        default = 0.14

        # matchcorr computation
        "--mxrot"
        help = "Maximum rotation"
        arg_type = Int64
        default = 10

        "--psi"
        help = "Minimum psi-s correlation"
        arg_type = Float64
        default = 0.95

        "--sz"
        help = "Minimum side length of floe mask"
        arg_type = Int64
        default = 16

        "--comp"
        help = "Size comparability"
        arg_type = Float64
        default = 0.25

        "--mm"
        help = "Maximum registration mismatch"
        arg_type = Float64
        default = 0.22

        # Goodness of match
        "--corr"
        help = "Mininimun psi-s correlation"
        arg_type = Float64
        default = 0.68

        "--area2"
        help = "Area thresholds to use for pairing floes"
        arg_type = Float64
        default = 0.236

        "--area3"
        help = "Area thresholds to use for pairing floes"
        arg_type = Float64
        default = 0.18
    end
    return nothing
end

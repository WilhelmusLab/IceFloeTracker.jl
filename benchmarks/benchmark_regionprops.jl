#=
IceFloeTracker `regionprops` vs scikit-image `regionprops_table` -- JULIA SIDE.

Run:
    julia --project=benchmarks -t 1 benchmarks/benchmark_regionprops.jl
    julia --project=benchmarks -t 1 benchmarks/benchmark_regionprops.jl --samples 3

Requires the fixtures from `export_labels.jl`. Emits
`benchmarks/results/julia.json` (timings) and
`benchmarks/results/values_julia.csv` (computed property values, for the parity
join). One scene, every floe in it -- no size sweep. Nothing here compares against Python; `join_results.jl` does that.

WHAT IS MEASURED

Only the properties that map 1:1 onto a scikit-image property:

    area, perimeter, convex_area, solidity, major_axis_length, minor_axis_length

`PolygonConvexArea` is excluded from this benchmark entirely. It returns a
continuous polygon area over pixel centres, which is not commensurable with a
pixel count: `:solidity` built on it exceeds 1 for 57% of real floes, and
scikit-image has no polygon-integration variant to compare it against. Timing a
quantity that is both defective and unmatched would add nothing. See
`issue_polygon_solidity.md`. Every measurement here pins
`convex_area_algorithm=PixelConvexArea()`.

`circularity` is excluded because it has no scikit-image counterpart (it is
`area/perimeter` here, not the dimensionless 4piA/P^2). `centroid`/`bbox` differ
by indexing convention only, `orientation` by sign and reference axis, and
`mask` is an array rather than a scalar -- all four are out of scope for this
study and are listed as such in the report.

HOW THE COMPARISON IS KEPT FAIR

  * Both sides read the SAME Int32 array from disk and assert its SHA-256, so a
    difference can never be an input difference.
  * Each property is timed in isolation AND the whole set is timed in one call.
    Both numbers are needed: `regionprops` shares `component_lengths`,
    `component_boxes` and the moment pass across properties, so the per-property
    timings do NOT sum to the combined one. Reporting only the sum would
    overstate the cost of a realistic call; reporting only the combined figure
    would hide which property dominates.
  * Every measured call is warmed up once first, so no compilation lands in a
    sample.
  * Inputs are interpolated with `$` into the benchmark expression. Without it
    the measurement includes a global-variable lookup rather than just the call.
  * `minimum` is the reported estimator, and the Python side is read the same
    way. Comparing Julia's minimum against pytest-benchmark's mean would flatter
    Julia for free.
  * GC is left enabled on both sides (`gcsample=false` here, and
    pytest-benchmark does not disable it either). Allocation figures are Julia
    only -- there is no symmetric Python number, so they get their own column
    rather than being passed off as comparable.
  * `-t 1` and `OMP_NUM_THREADS=1` on both sides.

ALIGNMENT TRAP THIS SCRIPT ENCODES

`regionprops` keeps labels with `area > minimum_area` -- strictly greater,
default 1 -- so single-pixel regions are dropped. scikit-image keeps them. The
value CSV written here therefore contains only the labels Julia kept, and the
Python side applies the identical filter. If the two row counts ever disagree,
the join refuses to run rather than silently comparing misaligned vectors.
=#

using BenchmarkTools
using DataFrames: nrow
using Printf: @printf
using SHA: sha256
using IceFloeTracker: regionprops_table, PixelConvexArea
import ImageMorphology

const BENCH_DIR = @__DIR__
const FIXTURE_DIR = joinpath(BENCH_DIR, "fixtures")
const RESULTS_DIR = joinpath(BENCH_DIR, "results")

# The properties that have a like-for-like scikit-image counterpart. Order is
# fixed so the combined-call case is reproducible.
const PROPERTIES = [
    :area, :perimeter, :convex_area, :solidity, :major_axis_length, :minor_axis_length
]

# Every property `regionprops_table` computes by default, plus the two derived
# ones. Used by `--properties all`, which is for measuring a change inside
# IceFloeTracker or its dependencies rather than comparing against scikit-image;
# `:mask`, `:bbox`, `:centroid` and `:orientation` have no like-for-like
# scikit-image counterpart and are excluded from PROPERTIES for that reason.
#
# `PROPERTIES_REQUIRING_BBOXES` in regionprops.jl is {bbox, perimeter,
# convex_area, mask}; solidity and circularity inherit the requirement through
# convex_area and perimeter. Those six are the ones a `component_boxes` change
# can move, and the rest act as controls.
const ALL_PROPERTIES = [
    :area, :perimeter, :convex_area, :solidity, :circularity, :major_axis_length,
    :minor_axis_length, :centroid, :bbox, :orientation, :mask,
]

"""
    apply_patch(path)

Evaluate the `component_boxes` definition from `path` (a checkout of
ImageMorphology's `src/connected.jl`) into the loaded ImageMorphology, replacing
that one method and nothing else.

This is how the pending JuliaImages/ImageMorphology.jl#146 is measured
end-to-end. `Pkg.develop` of that clone is not an option here: its master raises
the DataStructures bound beyond what this project's TiffImages pin allows, so it
cannot be resolved into this environment.
"""
function apply_patch(path)
    src = read(path, String)
    m = match(r"(?ms)^function component_boxes\(A::AbstractArray\{T,N\}.*?^end$", src)
    isnothing(m) && error("no component_boxes definition in $path")
    Base.eval(ImageMorphology, Meta.parse(m.match))
    println("patched ImageMorphology.component_boxes from ", path)
    return nothing
end

function argvalue(flag, default)
    i = findfirst(==(flag), ARGS)
    isnothing(i) && return default
    i == length(ARGS) && error("$flag requires a value")
    return ARGS[i + 1]
end

# ---------------------------------------------------------------- fixture I/O

function read_sidecar(path)
    d = Dict{String,Any}()
    for line in eachline(path)
        m = match(r"\"(\w+)\":\s*(.+?),?\s*$", line)
        isnothing(m) && continue
        k, raw = m.captures[1], m.captures[2]
        d[k] = startswith(raw, "\"") ? String(strip(raw, ['"'])) : parse(Int, raw)
    end
    return d
end

"""
Read a fixture written by `export_labels.jl`, verifying its digest.

The file is C-ordered (rows-major) so NumPy can read it without a transpose;
Julia is column-major, hence reading into a (cols, rows) buffer and permuting.
"""
function load_fixture(name)
    binpath = joinpath(FIXTURE_DIR, "labels_$name.bin")
    metapath = joinpath(FIXTURE_DIR, "labels_$name.json")
    isfile(binpath) || error("missing fixture $binpath -- run export_labels.jl first")
    meta = read_sidecar(metapath)

    digest = bytes2hex(open(sha256, binpath))
    digest == meta["sha256"] || error(
        "fixture $name digest mismatch:\n  on disk: $digest\n  sidecar: $(meta["sha256"])\n" *
        "regenerate with export_labels.jl",
    )

    rows, cols = meta["rows"], meta["cols"]
    buf = Vector{Int32}(undef, rows * cols)
    read!(binpath, buf)
    labels = permutedims(reshape(ltoh.(buf), cols, rows))
    return Matrix{Int64}(labels), meta
end

# ---------------------------------------------------------------- reporting

const RESULTS = Vector{Dict{String,Any}}()
const TAG = Ref("")

function record(case, property, meta, trial)
    stats = Dict{String,Any}(
        "case" => case,
        "impl" => "julia",
        "tag" => TAG[],
        "property" => property,
        "size" => meta["name"],
        "n_px" => meta["rows"] * meta["cols"],
        "n_labels" => meta["n_labels"],
        "min_s" => minimum(trial).time / 1e9,
        "median_s" => BenchmarkTools.median(trial).time / 1e9,
        "median_mib" => BenchmarkTools.median(trial).memory / 1024^2,
        "samples" => length(trial.times),
    )
    push!(RESULTS, stats)
    @printf("  %-34s min %8.4f s   median %8.4f s   %8.1f MiB  (n=%d)\n",
            property, stats["min_s"], stats["median_s"], stats["median_mib"],
            stats["samples"])
    return stats
end

function git_describe()
    try
        commit = readchomp(`git -C $(dirname(BENCH_DIR)) rev-parse --short HEAD`)
        branch = readchomp(`git -C $(dirname(BENCH_DIR)) rev-parse --abbrev-ref HEAD`)
        return "$branch@$commit"
    catch
        return "unknown"
    end
end

function write_json(path, rows, env)
    open(path, "w") do io
        println(io, "{")
        println(io, "  \"env\": {")
        ks = sort(collect(keys(env)))
        for (i, k) in enumerate(ks)
            println(io, "    \"$k\": \"$(env[k])\"", i == length(ks) ? "" : ",")
        end
        println(io, "  },")
        println(io, "  \"results\": [")
        for (i, r) in enumerate(rows)
            fields = String[]
            for k in sort(collect(keys(r)))
                v = r[k]
                push!(fields, "\"$k\": " * (v isa AbstractString ? "\"$v\"" : string(v)))
            end
            println(io, "    {", join(fields, ", "), "}", i == length(rows) ? "" : ",")
        end
        println(io, "  ]")
        println(io, "}")
    end
    return path
end

"""
Dump the property values Julia computed, one row per kept label.

`solidity` and `convex_area` are NaN for regions below the convex-area minimum
(4 px) and where the hull degenerates. They are written as-is; the join counts
them as non-comparable rather than dropping or averaging them.
"""
function write_values(path, labels)
    # `convex_area_algorithm` is passed explicitly rather than left to the default
    # so that a future change of default cannot silently pull a polygon-derived
    # `convex_area`/`solidity` into the parity table.
    df = regionprops_table(
        labels;
        properties=vcat(:label, PROPERTIES),
        minimum_area=1,
        convex_area_algorithm=PixelConvexArea(),
    )
    cols = vcat(:label, PROPERTIES)
    open(path, "w") do io
        println(io, join(string.(cols), ","))
        for i in 1:nrow(df)
            println(io, join([string(df[i, c]) for c in cols], ","))
        end
    end
    return nrow(df)
end

# ---------------------------------------------------------------- main

function main()
    samples = parse(Int, argvalue("--samples", "5"))
    seconds = parse(Float64, argvalue("--seconds", "600"))
    tag = argvalue("--tag", "")
    TAG[] = tag
    suffix = isempty(tag) ? "" : "_" * tag
    patch = argvalue("--patch", nothing)
    isnothing(patch) || apply_patch(patch)
    properties = argvalue("--properties", "skimage") == "all" ? ALL_PROPERTIES : PROPERTIES
    mkpath(RESULTS_DIR)

    println("IceFloeTracker regionprops_table -- Julia side")
    println("julia    ", VERSION)
    println("threads  ", Threads.nthreads())
    println("commit   ", git_describe())
    println()

    labels, meta = load_fixture("full")
    @printf("%d x %d (%.1f Mpx), %d labels\n",
            meta["rows"], meta["cols"], meta["rows"] * meta["cols"] / 1e6,
            meta["n_labels"])

    # `convex_area_algorithm` is pinned rather than left to the default, so
    # `:solidity` is always the pixel-counted ratio and cannot become the
    # polygon-derived one if the default ever moves.
    alg = PixelConvexArea()

    # Per-property, in isolation.
    for p in properties
        props = [p]
        regionprops_table(labels; properties=props, convex_area_algorithm=alg)  # warmup
        trial = @benchmark regionprops_table($labels; properties=$props,
                                             convex_area_algorithm=$alg) evals = 1 samples =
            samples seconds = seconds
        record("single", string(p), meta, trial)
    end

    # The whole set in one call -- what a real caller does. Shared passes are
    # amortised here, so this is NOT the sum of the rows above.
    allprops = properties
    regionprops_table(labels; properties=allprops, convex_area_algorithm=alg)
    trial = @benchmark regionprops_table($labels; properties=$allprops,
                                         convex_area_algorithm=$alg) evals = 1 samples =
        samples seconds = seconds
    record("combined", "all_properties", meta, trial)

    vpath = joinpath(RESULTS_DIR, "values_julia$suffix.csv")
    n = write_values(vpath, labels)
    @printf("  values: %d rows -> %s\n\n", n, basename(vpath))

    env = Dict(
        "julia" => string(VERSION),
        "threads" => string(Threads.nthreads()),
        "commit" => git_describe(),
        "estimator" => "minimum",
        "tag" => tag,
        "image_morphology" => string(pkgversion(ImageMorphology)),
        "patch" => isnothing(patch) ? "none" : patch,
    )
    out = write_json(joinpath(RESULTS_DIR, "julia$suffix.json"), RESULTS, env)
    println("wrote ", out, " (", length(RESULTS), " rows)")
end

main()

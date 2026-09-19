# Perimeter benchmark, isolating the `component_boxes` bottleneck it sits behind.
#
# `regionprops_table(:perimeter)` needs a bounding box per label, so it calls
# `ImageMorphology.component_boxes`. Until JuliaImages/ImageMorphology.jl#146
# lands, that function boxes a captured variable (JuliaLang/julia#15276) and
# dominates the call, hiding everything else. `--patch` measures the fixed
# version without waiting for a release; see `--patch` below.
#
# Reports, in order: `component_boxes` alone, `component_floes` (the per-floe
# crops), `component_perimeters` (the estimator), and the whole
# `regionprops_table(:perimeter)` call, so a change can be attributed to a stage.
#
#   julia --project=benchmarks benchmarks/bench_perimeter.jl [options]
#
#     --samples N     BenchmarkTools samples, default 5
#     --tag NAME      suffix for results/perimeter_NAME.{json,csv}, default "untagged"
#     --patch PATH    path to a checkout of ImageMorphology's src/connected.jl;
#                     its `component_boxes` is evaluated into the loaded module
#                     before timing. Used instead of `Pkg.develop` because the
#                     current ImageMorphology master raises its DataStructures
#                     bound beyond what this project's TiffImages pin allows, so
#                     the clone cannot be resolved into this environment. One
#                     method is replaced and nothing else.
#
using BenchmarkTools, ImageMorphology, IceFloeTracker, Printf, SHA, JSON

const BENCH_DIR = @__DIR__
const FIXTURE_DIR = joinpath(BENCH_DIR, "fixtures")
const RESULTS_DIR = joinpath(BENCH_DIR, "results")

function argvalue(flag, default)
    i = findfirst(==(flag), ARGS)
    isnothing(i) && return default
    return ARGS[i + 1]
end

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

function load_fixture(name)
    binpath = joinpath(FIXTURE_DIR, "labels_$name.bin")
    meta = read_sidecar(joinpath(FIXTURE_DIR, "labels_$name.json"))
    digest = bytes2hex(open(sha256, binpath))
    digest == meta["sha256"] || error("fixture digest mismatch")
    rows, cols = meta["rows"], meta["cols"]
    buf = Vector{Int32}(undef, rows * cols)
    read!(binpath, buf)
    return Matrix{Int64}(permutedims(reshape(ltoh.(buf), cols, rows))), meta
end

samples = parse(Int, argvalue("--samples", "5"))
tag = argvalue("--tag", "untagged")

# `--patch <path>` replaces `component_boxes` in the loaded ImageMorphology with
# the definition from a local checkout, so the fix can be measured end-to-end
# without a Pkg resolve (the clone's master bumps DataStructures beyond what
# this project's TiffImages pin allows). The pre-fix source in that checkout is
# byte-identical to the installed release, so this substitutes one method and
# nothing else.
patch = argvalue("--patch", nothing)
if patch !== nothing
    src = read(patch, String)
    m = match(r"(?ms)^function component_boxes\(A::AbstractArray\{T,N\}.*?^end$", src)
    isnothing(m) && error("no component_boxes definition in $patch")
    Base.eval(ImageMorphology, Meta.parse(m.match))
    println("patched ImageMorphology.component_boxes from ", patch)
end

labels, meta = load_fixture("full")
im_version = pkgversion(ImageMorphology)
im_path = dirname(dirname(pathof(ImageMorphology)))

@printf("tag              %s\n", tag)
@printf("ImageMorphology  %s  (%s)\n", im_version, im_path)
@printf("julia            %s\n", VERSION)
@printf("scene            %d x %d, %d labels\n\n",
        meta["rows"], meta["cols"], meta["n_labels"])

rows = Dict{String,Any}[]

function bench(name, f)
    f()  # warmup
    t = @benchmark $f() evals = 1 samples = samples seconds = 900
    r = Dict{String,Any}(
        "tag" => tag, "case" => name,
        "image_morphology" => string(im_version),
        "min_s" => minimum(t).time / 1e9,
        "median_s" => BenchmarkTools.median(t).time / 1e9,
        "median_mib" => BenchmarkTools.median(t).memory / 1024^2,
        "samples" => length(t.times),
    )
    push!(rows, r)
    @printf("  %-28s min %8.4f s   median %8.4f s   %9.1f MiB  (n=%d)\n",
            name, r["min_s"], r["median_s"], r["median_mib"], r["samples"])
    return r
end

bench("component_boxes", () -> ImageMorphology.component_boxes(labels))
bench("component_floes", () -> IceFloeTracker.component_floes(labels))
bench("component_perimeters", () -> IceFloeTracker.component_perimeters(labels))
bench("regionprops_perimeter",
      () -> regionprops_table(labels; properties=[:perimeter], minimum_area=1))

# Value parity: the fix must not change a single perimeter.
df = regionprops_table(labels; properties=[:label, :perimeter], minimum_area=1)
vpath = joinpath(RESULTS_DIR, "perimeter_values_$tag.csv")
open(vpath, "w") do io
    println(io, "label,perimeter")
    for i in 1:size(df, 1)
        println(io, df[i, :label], ",", df[i, :perimeter])
    end
end
@printf("\n  values: %d rows -> %s\n", size(df, 1), basename(vpath))

mkpath(RESULTS_DIR)
open(joinpath(RESULTS_DIR, "perimeter_$tag.json"), "w") do io
    JSON.print(io, Dict("env" => Dict("julia" => string(VERSION),
                                      "image_morphology" => string(im_version),
                                      "image_morphology_path" => im_path,
                                      "tag" => tag),
                        "results" => rows), 2)
end
println("wrote results/perimeter_$tag.json")

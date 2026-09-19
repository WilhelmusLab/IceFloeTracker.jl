#=
Join the Julia and Python regionprops benchmark outputs into one markdown
report, ready to paste into a GitHub Discussion.

Run:
    julia --project=benchmarks -t 1 benchmarks/join_results.jl
    julia --project=benchmarks -t 1 benchmarks/join_results.jl --out report.md

Reads `results/julia.json`, `results/python.json` (pytest-benchmark format) and
the per-size `values_{julia,python}_<size>.csv` pairs. Writes
`results/report.md`.

TWO RULES THIS SCRIPT ENFORCES

 1. Same estimator on both sides. It reads pytest-benchmark's `min`, never its
    `mean`, because the Julia side reports `minimum`. Mixing the two would
    flatter Julia by roughly the width of its noise distribution.
 2. No speed number without a parity verdict next to it. A property that
    computes something different is not faster, it is different, and the tables
    are ordered parity-first for that reason.
=#

using JSON
using Printf: @printf, @sprintf

const BENCH_DIR = @__DIR__
const RESULTS_DIR = joinpath(BENCH_DIR, "results")

const PROPERTIES = [
    "area", "perimeter", "convex_area", "solidity", "major_axis_length", "minor_axis_length"
]

# Per-property tolerance for the parity verdict. `area` is an integer pixel
# count and must agree exactly; the rest are floating-point reductions whose
# operation order legitimately differs between implementations, so they get a
# relative tolerance at roughly the level double-precision accumulation
# explains. Anything above this is an algorithmic difference, not arithmetic.
const RTOL = Dict(
    "area" => 0.0,
    "perimeter" => 1e-9,
    "convex_area" => 1e-9,
    "solidity" => 1e-9,
    "major_axis_length" => 1e-9,
    "minor_axis_length" => 1e-9,
)

const EXCLUSIONS = [
    ("circularity", "No scikit-image counterpart. Implemented here as `area/perimeter`, which is not the dimensionless `4piA/P^2` either, so there is nothing to compare it against."),
    ("centroid", "Differs by indexing convention only (1-based inclusive vs 0-based). Comparable in principle; out of scope for this study."),
    ("bbox", "As `centroid`: 1-based inclusive vs 0-based half-open."),
    ("orientation", "Sign and reference axis differ (`0.5*atan(2mu11, mu20-mu02)` vs scikit-image's convention). Reconciling it is its own study."),
    ("mask", "A cropped `BitMatrix` vs scikit-image's `image`; an array, not a scalar to diff."),
    ("convex_area via PolygonConvexArea", "A continuous polygon area over pixel centres, not commensurable with a pixel count (solidity built on it exceeds 1 for 57% of real floes), and scikit-image has no polygon-integration variant to compare against. All measurements use `PixelConvexArea`."),
]

function argvalue(flag, default)
    i = findfirst(==(flag), ARGS)
    isnothing(i) && return default
    i == length(ARGS) && error("$flag requires a value")
    return ARGS[i + 1]
end

# ---------------------------------------------------------------- loading

"""Julia timings keyed by (case, property, size)."""
function load_julia(path)
    doc = JSON.parsefile(path)
    out = Dict{Tuple{String,String,String},Any}()
    for r in doc["results"]
        out[(r["case"], r["property"], r["size"])] = r
    end
    return out, doc["env"]
end

"""
Python timings keyed the same way, from either harness.

Two formats are accepted, distinguished by their top-level key:

  * `results` -- the `timeit` notebook (`03-python-timeit.ipynb`), same shape as
    `julia.json`.
  * `benchmarks` -- `bench_regionprops.py` under pytest-benchmark, where the
    case/property/size triple travels in `extra_info`.

Either way the join is on an identical `(case, property, size)` key rather than
on parsed test-id strings. `min` is read in both cases, never `mean`: the Julia
side reports `minimum`, and mixing estimators would flatter Julia for free.
"""
function load_python(path)
    doc = JSON.parsefile(path)
    out = Dict{Tuple{String,String,String},Any}()

    if haskey(doc, "results")               # timeit notebook
        for r in doc["results"]
            out[(r["case"], r["property"], r["size"])] = Dict(
                "min_s" => r["min_s"],
                "median_s" => r["median_s"],
                "samples" => r["samples"],
                "n_labels" => r["n_labels"],
                "n_px" => r["n_px"],
            )
        end
        env = get(doc, "env", Dict())
        return out, Dict(
            "python_version" => get(env, "python", "?"),
            "cpu" => Dict("brand_raw" => get(env, "machine", "?")),
            "timer" => get(env, "timer", "timeit"),
        )
    end

    for b in doc["benchmarks"]              # pytest-benchmark
        info = b["extra_info"]
        haskey(info, "case") || continue
        out[(info["case"], info["property"], info["size"])] = Dict(
            "min_s" => b["stats"]["min"],
            "median_s" => b["stats"]["median"],
            "samples" => b["stats"]["rounds"],
            "n_labels" => info["n_labels"],
            "n_px" => info["n_px"],
        )
    end
    machine = get(doc, "machine_info", Dict())
    machine["timer"] = "pytest-benchmark"
    return out, machine
end

function read_values(path)
    isfile(path) || return nothing
    lines = readlines(path)
    header = split(lines[1], ",")
    cols = Dict(h => Float64[] for h in header)
    labels = Int[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        parts = split(line, ",")
        push!(labels, parse(Int, parts[1]))
        for (h, p) in zip(header, parts)
            h == "label" && continue
            push!(cols[h], parse(Float64, p))
        end
    end
    return labels, cols
end

"""
Render property names as an English list of code spans: "`a`, `b` and `c`".
"""
function proplist(props)
    spans = ["`$p`" for p in props]
    length(spans) == 1 && return spans[1]
    return join(spans[1:(end - 1)], ", ") * " and " * spans[end]
end

# ---------------------------------------------------------------- parity

"""
Compare one property's values across the two implementations.

Rows where either side is non-finite are counted as non-comparable rather than
dropped silently or folded into the statistics: Julia's convex-area returns NaN
below 4 px and on hull degeneracy, where scikit-image returns a number, and
that gap is a finding in its own right.
"""
function compare_property(jv, pv, prop)
    n = length(jv)
    noncomparable = 0
    maxabs = 0.0
    maxrel = 0.0
    worst_label = 0
    compared = 0
    for i in 1:n
        a, b = jv[i], pv[i]
        if !isfinite(a) || !isfinite(b)
            noncomparable += 1
            continue
        end
        compared += 1
        d = abs(a - b)
        scale = max(abs(a), abs(b))
        r = scale == 0 ? 0.0 : d / scale
        if r > maxrel
            maxrel = r
            worst_label = i
        end
        maxabs = max(maxabs, d)
    end
    tol = RTOL[prop]
    pass = compared > 0 && maxrel <= tol
    return (; compared, noncomparable, maxabs, maxrel, pass, worst_label)
end

fmt_s(x) = x >= 1 ? @sprintf("%.3f s", x) :
           x >= 1e-3 ? @sprintf("%.2f ms", x * 1e3) : @sprintf("%.1f us", x * 1e6)

function main()
    outpath = argvalue("--out", joinpath(RESULTS_DIR, "report.md"))
    jpath = joinpath(RESULTS_DIR, "julia.json")
    # The timeit notebook's output wins when both are present; it is the
    # pytest-free harness and the one kept current.
    tpath = joinpath(RESULTS_DIR, "python_timeit.json")
    ppath = isfile(tpath) ? tpath : joinpath(RESULTS_DIR, "python.json")
    isfile(jpath) || error("missing $jpath -- run benchmark_regionprops.jl first")
    isfile(ppath) || error(
        "missing a Python result file -- run 03-python-timeit.ipynb (writes " *
        "python_timeit.json) or bench_regionprops.py (writes python.json)")
    println("python timings: ", basename(ppath))

    jul, jenv = load_julia(jpath)
    pyt, machine = load_python(ppath)

    io = IOBuffer()
    println(io, "# `regionprops_table`: IceFloeTracker.jl vs scikit-image\n")
    println(io, "Generated by `benchmarks/join_results.jl`. ",
                "Reproduce with `export_labels.jl`, `benchmark_regionprops.jl`, ",
                "then `bench_regionprops.py`.\n")

    # ---- environment
    println(io, "## Environment\n")
    println(io, "| | |")
    println(io, "|---|---|")
    println(io, "| Julia | $(get(jenv, "julia", "?")) |")
    println(io, "| IceFloeTracker commit | `$(get(jenv, "commit", "?"))` |")
    if haskey(machine, "python_version")
        println(io, "| Python | $(machine["python_version"]) |")
    end
    if haskey(machine, "cpu") && haskey(machine["cpu"], "brand_raw")
        println(io, "| CPU | $(machine["cpu"]["brand_raw"]) |")
    end
    println(io, "| Threads | 1 (`julia -t 1`, `OMP_NUM_THREADS=1`) |")
    println(io, "| Estimator | minimum, both sides |")
    haskey(machine, "timer") && println(io, "| Python timer | $(machine["timer"]) |")
    println(io)
    println(io, "One scene, every floe in it: ",
                "`001-fram_strait-20120412.aqua.labeled.png`, 5680x3392, 1910 labels. ",
                "Both languages call `regionprops_table` on one shared, ",
                "SHA-256-verified `int32` fixture.\n")

    # ---- parity, computed before the findings so those can be derived from it
    parity = Dict{String,Any}()
    nfloes = 0
    jvals = read_values(joinpath(RESULTS_DIR, "values_julia.csv"))
    pvals = read_values(joinpath(RESULTS_DIR, "values_python.csv"))
    have_values = !(isnothing(jvals) || isnothing(pvals))
    if have_values
        jlabels, jcols = jvals
        plabels, pcols = pvals
        jlabels == plabels || error(
            "label sets differ ($(length(jlabels)) Julia vs $(length(plabels)) Python) " *
            "-- comparison refused")
        nfloes = length(jlabels)
        for prop in PROPERTIES
            parity[prop] = compare_property(jcols[prop], pcols[prop], prop)
        end
    end

    # ---- headline findings, computed from the tables below, never asserted
    # Anything here that is not derived from the measured rows does not belong here.
    ratios = Dict{String,Vector{Float64}}()
    for ((case, prop, size), j) in jul
        case == "algorithm" && continue
        haskey(pyt, (case, prop, size)) || continue
        push!(get!(ratios, prop, Float64[]), pyt[(case, prop, size)]["min_s"] / j["min_s"])
    end
    if !isempty(ratios)
        allr = collect(Iterators.flatten(values(ratios)))
        faster = count(>(1.0), allr)
        worst = argmin(p -> minimum(ratios[p]), collect(keys(ratios)))
        best = argmax(p -> maximum(ratios[p]), collect(keys(ratios)))
        maxalloc = maximum(r["median_mib"] for r in values(jul))

        println(io, "## Findings\n")
        println(io, "Every number here is computed from the tables below.\n")
        println(io, "- scikit-image is faster for **$(length(allr) - faster) of $(length(allr))** ",
                    "measured calls; Julia is faster for **$faster**.")
        println(io, "- Slowest relative to scikit-image: `$worst` at ",
                    @sprintf("%.2fx", minimum(ratios[worst])),
                    ". Closest: `$best` at ", @sprintf("%.2fx", maximum(ratios[best])), ".")
        matching = [p for p in PROPERTIES if haskey(parity, p) && parity[p].pass]
        differing = [p for p in PROPERTIES if haskey(parity, p) && !parity[p].pass]
        if !isempty(matching)
            println(io, "- ", proplist(matching), " agree with scikit-image to machine ",
                        "precision (max relative difference ",
                        @sprintf("%.0e", maximum(parity[p].maxrel for p in matching)), ").")
        end
        if !isempty(differing)
            println(io, "- ", proplist(differing), " **do not** agree, by up to ",
                        @sprintf("%.0f%%", 100 * maximum(parity[p].maxrel for p in differing)),
                        ". See `benchmarks/convex_area_discrepancy.md` for a minimal ",
                        "reproducer.")
        end
        behind = sort([p for p in keys(ratios) if minimum(ratios[p]) < 1.0];
                      by=p -> minimum(ratios[p]))
        if isempty(behind)
            println(io, "- Julia is at least as fast as scikit-image on every measured call.")
        else
            slowest = minimum(minimum(ratios[p]) for p in behind)
            agree = all(haskey(parity, p) && parity[p].pass for p in behind)
            tail = if agree
                ". Those values match, so what is left there is implementation, not method."
            else
                "."
            end
            println(io, "- scikit-image is still ahead on ", proplist(behind), ", down to ",
                        @sprintf("%.2fx", slowest), tail)
        end
        println(io, "- Peak Julia allocation in this sweep: ",
                    @sprintf("%.1f", maxalloc), " MiB.")
        println(io)
    end

    # ---- parity first, deliberately
    println(io, "## Parity\n")
    println(io, "Read this before the timings. A property that computes something ",
                "different is not faster, it is different.\n")
    println(io, "Rows are the labels `regionprops_table` keeps (`area > 1`); the Python side ",
                "applies the identical filter. \"n/c\" counts rows where one side is ",
                "non-finite -- Julia's convex area is `NaN` below 4 px and on hull ",
                "degeneracy, where scikit-image still returns a number.\n")

    if !have_values
        println(io, "*Value files missing -- run both harnesses.*\n")
    else
        println(io, "$nfloes floes compared.\n")
        println(io, "| property | compared | n/c | max abs diff | max rel diff | verdict |")
        println(io, "|---|---:|---:|---:|---:|:---:|")
        for prop in PROPERTIES
            r = parity[prop]
            println(io, "| `$prop` | $(r.compared) | $(r.noncomparable) | ",
                        @sprintf("%.3g", r.maxabs), " | ",
                        @sprintf("%.3g", r.maxrel), " | ",
                        r.pass ? "match" : "**differs**", " |")
        end
        println(io)
    end

    # ---- speed
    println(io, "## Speed\n")
    println(io, "Minimum over the samples, both sides. `ratio` is Python / Julia: ",
                "greater than 1 means Julia is faster.\n")
    println(io, "Per-property rows time that property **alone**. They do not sum to ",
                "`all_properties`: `regionprops` shares the label-lengths, bounding-box ",
                "and moment passes across properties, so a combined call amortises work ",
                "that each isolated call repeats.\n")

    rows = String[]
    for (case, prop) in vcat([("single", p) for p in PROPERTIES],
                             [("combined", "all_properties")])
        k = (case, prop, "full")
        (haskey(jul, k) && haskey(pyt, k)) || continue
        j, p = jul[k]["min_s"], pyt[k]["min_s"]
        note = case == "single" && haskey(parity, prop) && !parity[prop].pass ? " [^d]" : ""
        push!(rows, "| `$prop`$note | $(fmt_s(j)) | $(fmt_s(p)) | " *
                    @sprintf("%.2fx", p / j) * " | " *
                    @sprintf("%.1f", jul[k]["median_mib"]) * " |")
    end
    if !isempty(rows)
        println(io, "| property | Julia (min) | scikit-image (min) | ratio | Julia alloc (MiB) |")
        println(io, "|---|---:|---:|---:|---:|")
        println(io, join(rows, "\n"))
        println(io)
    end
    println(io, "[^d]: parity differs for this property -- see above; the ratio compares ",
                "two different computations.\n")

    # ---- exclusions
    println(io, "## Properties excluded from this study\n")
    println(io, "| property | why |")
    println(io, "|---|---|")
    for (p, why) in EXCLUSIONS
        println(io, "| `$p` | $why |")
    end
    println(io)

    # ---- method
    println(io, "## Method and caveats\n")
    println(io, """
- **One shared input.** Labels are decoded once by `export_labels.jl`, written as raw little-endian `int32`, and SHA-256'd. Both harnesses assert the digest before measuring, so a difference can never be an input difference.
- **Same estimator.** Julia's `minimum` against pytest-benchmark's `min`. Comparing it to pytest-benchmark's *mean* would have flattered Julia for free.
- **Warmup.** Both sides call once before timing; Julia's first call includes compilation and is excluded. Time-to-first-call is a separate question and is not measured here.
- **GC** is enabled on both sides (`gcsample=false`; pytest-benchmark does not disable it, unlike `timeit`).
- **`PolygonConvexArea` is excluded entirely.** Every measurement pins `convex_area_algorithm=PixelConvexArea()`. The polygon variant returns a continuous area over pixel centres, not commensurable with a pixel count: `:solidity` built on it exceeds 1 for 57% of real floes, and scikit-image has no polygon-integration variant to compare it to.
- **Allocations are Julia-only.** There is no symmetric scikit-image figure, so the column is labelled as such rather than presented as a comparison.
- **Threads** pinned to 1 on both sides, same machine, same session.
- **Row filtering.** `regionprops_table` keeps `area > minimum_area` with `minimum_area=1` by default, i.e. *strictly* greater, dropping single-pixel regions; scikit-image keeps everything. The Python harness applies the identical filter, and the join refuses to run if the label sets disagree.
- **One scene, all of it.** There is no size sweep: cropping changes the floe population as much as the pixel count, so a crop measures a different workload rather than a smaller one.
""")

    mkpath(RESULTS_DIR)
    write(outpath, String(take!(io)))
    println("wrote ", outpath)

    # console summary
    println("\nparity summary:")
    for prop in PROPERTIES
        haskey(parity, prop) || continue
        r = parity[prop]
        @printf("  %-18s %s  (maxrel %.3g, n/c %d)\n",
                prop, r.pass ? "match " : "DIFFER", r.maxrel, r.noncomparable)
    end
end

main()

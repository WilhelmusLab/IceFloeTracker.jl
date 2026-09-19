#=
Export the shared label-array fixture used by BOTH the Julia and the Python
regionprops benchmarks.

Run:
    julia --project=benchmarks -t 1 benchmarks/export_labels.jl

WHY A FILE INSTEAD OF EACH SIDE LOADING THE PNG

The comparison is only meaningful if both implementations are handed the *same*
integers. Two independent PNG loads (Images.jl vs imageio/PIL) agreeing is an
assumption, not a fact: 16-bit grayscale PNGs have gone through
fixed-point/float round-trips in either library before now. So the labels are
decoded once, here, written as raw little-endian Int32, and SHA-256'd. Both
harnesses assert the digest before they measure anything; a stale or
re-generated fixture then fails loudly instead of producing a silently
mismatched comparison.

Raw binary + a JSON sidecar keeps both readers dependency-free: `np.fromfile`
on one end, `read!` on the other. No NPZ/HDF5 dependency is added to either
environment just to move a rectangle of integers.

The whole scene is exported, with every floe in it. There is no size sweep:
cropping changes the floe population as much as the pixel count (the corner of
this scene is sparser than its middle, and a third of the floes in a small crop
are cut by the crop edge), so a crop measures a different workload rather than
a smaller one.
=#

using Images: load, gray
using SHA: sha256
using Printf: @printf

const BENCH_DIR = @__DIR__
const FIXTURE_DIR = joinpath(BENCH_DIR, "fixtures")
const DEFAULT_SCENE = joinpath(BENCH_DIR, "001-fram_strait-20120412.aqua.labeled.png")

function argvalue(flag, default)
    i = findfirst(==(flag), ARGS)
    isnothing(i) && return default
    i == length(ARGS) && error("$flag requires a value")
    return ARGS[i + 1]
end

"""Decode the 16-bit grayscale labeled PNG into an Int32 label matrix."""
function load_labels(path)
    isfile(path) || error("no such scene: $path")
    img = load(path)
    return Int32.(reinterpret.(UInt16, gray.(img)))
end

"""Number of distinct non-background labels present."""
n_labels(labels) = count(!=(0), unique(labels))

function json_escape(s)
    return replace(string(s), "\\" => "\\\\", "\"" => "\\\"")
end

function write_sidecar(path, d)
    open(path, "w") do io
        println(io, "{")
        keys_sorted = sort(collect(keys(d)))
        for (i, k) in enumerate(keys_sorted)
            v = d[k]
            rendered = v isa AbstractString ? "\"$(json_escape(v))\"" : string(v)
            print(io, "  \"$k\": $rendered")
            println(io, i == length(keys_sorted) ? "" : ",")
        end
        println(io, "}")
    end
    return path
end

"""
Write `labels` as raw little-endian Int32 plus a JSON sidecar, and return the
sidecar's metadata dict.

Stored row-major (C order) so NumPy can `fromfile(...).reshape(rows, cols)`
directly. Julia is column-major, hence the transpose on write and on read --
doing it here, once, keeps the orientation question out of both harnesses.
"""
function export_case(name, labels, source)
    mkpath(FIXTURE_DIR)
    binpath = joinpath(FIXTURE_DIR, "labels_$name.bin")
    rows, cols = size(labels)

    open(binpath, "w") do io
        write(io, htol.(permutedims(labels)))
    end

    digest = bytes2hex(open(sha256, binpath))
    meta = Dict(
        "name" => name,
        "rows" => rows,
        "cols" => cols,
        "dtype" => "int32",
        "byte_order" => "little",
        "order" => "C",
        "n_labels" => n_labels(labels),
        "sha256" => digest,
        "source" => basename(source),
    )
    write_sidecar(joinpath(FIXTURE_DIR, "labels_$name.json"), meta)

    @printf("%-8s %5d x %-5d  %7.1f Mpx  %5d labels  %s\n",
            name, rows, cols, rows * cols / 1e6, meta["n_labels"], digest[1:12])
    return meta
end

function main()
    scene = argvalue("--scene", DEFAULT_SCENE)

    println("scene: ", scene)
    labels = load_labels(scene)

    println("\nname         rows x cols        size   labels  sha256")
    export_case("full", labels, scene)

    println("\nfixture written to ", FIXTURE_DIR)
end

main()

# `regionprops` vs scikit-image — benchmark harness

Compares `IceFloeTracker.regionprops_table` against
`skimage.measure.regionprops_table` on a real Fram Strait scene — one scene,
all 1910 floes — for both **speed** and **numerical agreement**.

Results: `results/report.md` (paste-ready for a GitHub Discussion).

## Running it

Four steps, in order. Run them **sequentially on an otherwise idle machine** —
running the two harnesses at the same time makes them contend for CPU and
corrupts both sets of timings.

```bash
# 1. Build the shared fixture (once). Decodes the labeled PNG into raw int32
#    + a SHA-256 sidecar that both harnesses verify before measuring.
julia --project=benchmarks -t 1 benchmarks/export_labels.jl

# 2. Julia side -> results/julia.json + results/values_julia.csv
OMP_NUM_THREADS=1 julia --project=benchmarks -t 1 \
    benchmarks/benchmark_regionprops.jl --samples 10

# 3. Python side -> results/python.json + results/values_python.csv
cd benchmarks/python
OMP_NUM_THREADS=1 uv run pytest bench_regionprops.py \
    --benchmark-json=../results/python.json --benchmark-min-rounds=5
cd ../..

# 4. Join into the report
julia --project=benchmarks -t 1 benchmarks/join_results.jl
```

Both harnesses take well under a minute at `--samples 10`. Before the
optimisation the Julia side took roughly 12 minutes, nearly all of it the
convex-area and perimeter passes.

### Reproducing the before/after numbers

The committed results were produced against three states, because most of the
speedup depends on a change that is not in this repository:

| results file | state measured |
|---|---|
| `results/julia_allprops-baseline.json` | `main` before the optimisation |
| `results/julia_pr-nopatch.json` | this branch, stock ImageMorphology |
| `results/julia.json` | this branch + ImageMorphology.jl#146 |

[JuliaImages/ImageMorphology.jl#146](https://github.com/JuliaImages/ImageMorphology.jl/pull/146)
removes a captured-variable boxing in `component_boxes`
([JuliaLang/julia#15276](https://github.com/JuliaLang/julia/issues/15276)) that
otherwise dominates every bounding-box-dependent property. Until it is released,
`--patch` measures it without a `Pkg.develop` — that clone cannot be resolved
into this environment, because its master raises the DataStructures bound beyond
what this project's TiffImages pin allows:

```bash
git clone https://github.com/JuliaImages/ImageMorphology.jl ~/ImageMorphology.jl
# check out the PR branch, then:
julia --project=benchmarks -t 1 benchmarks/benchmark_regionprops.jl \
    --samples 10 --patch ~/ImageMorphology.jl/src/connected.jl
```

`apply_patch` evaluates that file's `component_boxes` into the loaded
ImageMorphology, replacing one method and nothing else.

### Other flags

`benchmark_regionprops.jl` also takes:

- `--tag NAME` — write `results/julia_NAME.json` and `results/values_julia_NAME.csv`
  instead of the default names, so a run does not overwrite the committed ones.
  `join_results.jl` always reads the untagged pair.
- `--properties all` — time every property `regionprops_table` computes by
  default plus `circularity`, rather than only those with a scikit-image
  counterpart. For attributing a change inside IceFloeTracker; the report is
  built from the default set.

`bench_perimeter.jl` breaks the perimeter call into `component_boxes`,
`component_floes`, `component_perimeters` and the whole `regionprops_table`
call, so a change can be attributed to a stage. It takes the same `--samples`,
`--tag` and `--patch`.

### Committed results

`results/` is ignored except for the files the pull request quotes: `report.md`,
the three `julia*.json` above, `python.json`, and both value CSVs. Runs from
developing the change are not kept. `fixtures/` is ignored too — step 1
reproduces it byte-for-byte, and the digest is checked on every read.

The Python environment is managed by `uv` and pinned by `python/uv.lock`;
`uv run` builds it on first use. The `.bin` fixtures are regenerable and are not
worth committing — step 1 reproduces them byte-for-byte (the digest is stable
across runs).

## Interactive exploration

Two notebooks in `notebooks/`, for poking at the benchmark without re-running
the whole sweep:

| notebook | kernel | covers |
|---|---|---|
| `01-julia-side.ipynb` | `julia-1.12` | load a fixture, browse the property values, re-time any single call, per-floe drill-down by area bucket |
| `02-comparison.ipynb` | project venv | the scikit-image side, the parity and speed joins, and a visual walk through the `convex_area` disagreement |
| `03-python-timeit.ipynb` | project venv | pytest-free Python timings using stdlib `timeit`; exports `results/python_timeit.json` |

```bash
# Python notebook, with the project's pinned scikit-image:
cd benchmarks/python && uv run jupyter lab ../notebooks/02-comparison.ipynb

# Julia notebook: any Jupyter with the julia-1.12 kernel. Its first cell
# activates benchmarks/Project.toml, so no manual environment setup.
jupyter lab benchmarks/notebooks/01-julia-side.ipynb
```

All three run on the full scene. They are committed with outputs so they read
as a report without being run.

Plots live only in the Python notebook — the `benchmarks` Julia environment has
no plotting dependency and this keeps it that way.

### If the Julia notebook fails to precompile

The symptom is several minutes of ordinary-looking precompile output ending in

```
ERROR: LoadError: UndefVarError: `StaticData` not defined in `Base`
... Failed to precompile BenchmarkTools ... /.julia/compiled/v1.11/...
```

`Base.StaticData` is Julia 1.12 only. IceFloeTracker requires `julia = "1.12"`
and `benchmarks/Manifest.toml` is resolved for it, so on a 1.11 kernel
`PrecompileTools` fails and takes `PrettyTables`, `DataFrames` and
`IceFloeTracker` down with it. The `v1.11` in the cache path is the tell.

**In VS Code this is the easy way to get it wrong**: the Julia extension
supplies its own notebook kernel and ignores the Jupyter kernelspec entirely,
so the notebook runs on whatever `julia.executablePath` points at. Set it to a
1.12 binary, or pick "Julia 1.12" in the kernel selector.

The notebook's first cell checks `VERSION` and stops immediately with these
instructions, rather than letting the cascade run.

Separately, if the Jupyter `julia-1.12` kernel fails to *start* after a
`juliaup` upgrade, its kernelspec may point at a version that no longer
exists. Reinstall it with `using IJulia; installkernel("Julia")`; pointing
`argv[0]` at `~/.juliaup/bin/julia` instead of a pinned path avoids the
recurrence.

## Files

| file | role |
|---|---|
| `export_labels.jl` | decodes the scene once, writes the shared `int32` fixture + SHA-256 sidecar |
| `benchmark_regionprops.jl` | Julia timings (BenchmarkTools) + computed values |
| `bench_perimeter.jl` | the perimeter call broken into its stages |
| `python/bench_regionprops.py` | scikit-image timings (pytest-benchmark) + computed values |
| `join_results.jl` | joins both, emits `results/report.md` |
| `convex_area_discrepancy.md` | the one parity failure, with a minimal reproducer |
| `notebooks/01-julia-side.ipynb` | interactive Julia-side exploration |
| `notebooks/02-comparison.ipynb` | interactive scikit-image side, parity, and the convex-area investigation |
| `notebooks/03-python-timeit.ipynb` | pytest-free Python timings with stdlib `timeit` |
| `issue_polygon_solidity{.md,_mwe.jl}` | separate `PolygonConvexArea` defect: `solidity > 1` |

## What is and isn't compared

Six properties have a like-for-like scikit-image counterpart and are studied:
`area`, `perimeter`, `convex_area`, `solidity`, `major_axis_length`,
`minor_axis_length`.

`PolygonConvexArea` is excluded too: it returns a continuous polygon area over
pixel centres, which is not commensurable with a pixel count, and scikit-image
has no polygon-integration variant to compare it against. Every measurement
pins `convex_area_algorithm=PixelConvexArea()`. See `issue_polygon_solidity.md`.

`circularity` (no counterpart — it is `area/perimeter` here, not `4piA/P^2`),
`centroid` and `bbox` (indexing convention only), `orientation` (sign and
reference axis) and `mask` (an array, not a scalar) are out of scope. The report
lists each with its reason.

## Two Python timing harnesses

Either can produce the timings; `join_results.jl` reads whichever is present and
prefers `python_timeit.json`, naming the timer in the report's Environment
table.

| | `bench_regionprops.py` | `03-python-timeit.ipynb` |
|---|---|---|
| tool | pytest-benchmark | stdlib `timeit` |
| output | `results/python.json` | `results/python_timeit.json` |
| rounds | auto-calibrated | `repeat=5, number=1` |
| GC | enabled (its default) | enabled **explicitly** — see below |

**`timeit` disables the garbage collector while timing.** That is its
documented behaviour and sensible for a microbenchmark, but BenchmarkTools
leaves GC on, so an uncorrected `timeit` number would flatter Python. The
notebook re-enables it from `setup`, which runs inside the timed callable and
therefore after `timeit`'s own `gc.disable()`:

```python
timeit.Timer(stmt, setup="import gc; gc.enable()")
```

Measured on `area_convex`, the GC correction is worth well under 1% here — but
it is free, and it stops the comparison depending on a default that happens not
to matter.

## Design notes

Things that will silently produce a wrong comparison if changed carelessly:

- **One shared input.** Both languages read the same SHA-256-verified file. Two
  independent PNG decodes agreeing is an assumption, not a fact.
- **`regionprops_table` keeps `area > minimum_area`** — *strictly* greater, default 1
  — so single-pixel regions are dropped; scikit-image keeps them. The Python
  harness applies the identical filter, and `join_results.jl` refuses to compare
  if the label sets disagree rather than lining up mismatched vectors. This
  filter removes 76 of 1910 labels, so it is not a hypothetical.
- **Same estimator on both sides**: Julia's `minimum` against pytest-benchmark's
  `min`, never its `mean`.
- **Per-property timings do not sum to `all_properties`.** `regionprops_table` shares
  the label-lengths, bounding-box and moment passes across properties.
- **Allocations are Julia-only** — there is no symmetric scikit-image figure, so
  the column is labelled rather than presented as a comparison.
- **One scene, all of it.** There is no size sweep. Cropping changes the floe
  population as much as the pixel count — the corner of this scene is sparser
  than its middle, and a third of the floes in a small crop are cut by the crop
  edge — so a crop measures a different workload rather than a smaller one.

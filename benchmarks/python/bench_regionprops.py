"""
IceFloeTracker `regionprops` vs scikit-image `regionprops_table` -- PYTHON SIDE.

Run (from benchmarks/python):
    OMP_NUM_THREADS=1 uv run pytest bench_regionprops.py \
        --benchmark-json=../results/python.json

One scene, every floe in it. Requires the fixture from `export_labels.jl`; the
SHA-256 in the sidecar is asserted before anything is measured, so this can
never time a different array from the one Julia timed.

Writes `../results/values_python.csv` alongside the timing JSON, for the parity
join.

The two alignment traps documented in benchmark_regionprops.jl are handled here:

  * `regionprops` keeps labels with `area > 1` (strictly greater); skimage keeps
    every label. `_filtered_table` applies the identical filter so the row sets
    match label-for-label.
  * Julia's convex-area returns NaN below 4 px and on hull degeneracy, whereas
    skimage always returns a number. Nothing is dropped here -- the join counts
    those rows as non-comparable.

Timing notes: pytest-benchmark calibrates its own rounds and does NOT disable
GC (unlike timeit), matching the Julia side's `gcsample=false`. The join reads
the `min` field from this JSON, not the mean, because the Julia side reports
`minimum`.
"""

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from skimage.measure import regionprops_table

BENCH_DIR = Path(__file__).resolve().parent.parent
FIXTURE_DIR = BENCH_DIR / "fixtures"
RESULTS_DIR = BENCH_DIR / "results"

# scikit-image names for the properties that map 1:1 onto the Julia ones.
# Keys are the Julia names so both sides emit identically-labelled cases.
PROPERTY_MAP = {
    "area": "area",
    "perimeter": "perimeter",
    "convex_area": "area_convex",
    "solidity": "solidity",
    "major_axis_length": "axis_major_length",
    "minor_axis_length": "axis_minor_length",
}

def load_fixture(name):
    """Read a fixture written by export_labels.jl, verifying its digest."""
    binpath = FIXTURE_DIR / f"labels_{name}.bin"
    metapath = FIXTURE_DIR / f"labels_{name}.json"
    if not binpath.exists():
        pytest.skip(f"missing fixture {binpath} -- run export_labels.jl first")
    meta = json.loads(metapath.read_text())

    digest = hashlib.sha256(binpath.read_bytes()).hexdigest()
    assert digest == meta["sha256"], (
        f"fixture {name} digest mismatch:\n  on disk: {digest}\n"
        f"  sidecar: {meta['sha256']}\n  regenerate with export_labels.jl"
    )

    # Written C-ordered little-endian int32 precisely so this is a plain read.
    labels = np.fromfile(binpath, dtype="<i4").reshape(meta["rows"], meta["cols"])
    return labels, meta


def _filtered_table(labels, properties):
    """
    regionprops_table restricted to the labels Julia would keep.

    `area` is always requested so the filter can be applied even when the
    caller only asked for, say, `perimeter`; it is dropped again afterwards
    unless it was asked for.
    """
    props = tuple(dict.fromkeys(("label", "area") + tuple(properties)))
    table = regionprops_table(labels, properties=props)
    keep = table["area"] > 1
    return {k: v[keep] for k, v in table.items()}


@pytest.mark.parametrize("julia_name", list(PROPERTY_MAP))
def test_single_property(benchmark, julia_name):
    """One property in isolation -- mirrors the Julia `single` case."""
    labels, meta = load_fixture("full")
    skname = PROPERTY_MAP[julia_name]

    benchmark.extra_info.update(
        {
            "case": "single",
            "impl": "python",
            "property": julia_name,
            "size": meta["name"],
            "n_px": meta["rows"] * meta["cols"],
            "n_labels": meta["n_labels"],
        }
    )
    benchmark(_filtered_table, labels, (skname,))


def test_all_properties(benchmark):
    """The whole set in one call -- mirrors the Julia `combined` case."""
    labels, meta = load_fixture("full")
    sknames = tuple(PROPERTY_MAP.values())

    benchmark.extra_info.update(
        {
            "case": "combined",
            "impl": "python",
            "property": "all_properties",
            "size": meta["name"],
            "n_px": meta["rows"] * meta["cols"],
            "n_labels": meta["n_labels"],
        }
    )
    result = benchmark(_filtered_table, labels, sknames)

    # Dump the values once, outside the timed region, for the parity join.
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out = RESULTS_DIR / "values_python.csv"
    cols = ["label"] + list(PROPERTY_MAP)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        n = len(result["label"])
        for i in range(n):
            w.writerow(
                [result["label"][i]]
                + [result[PROPERTY_MAP[p]][i] for p in PROPERTY_MAP]
            )

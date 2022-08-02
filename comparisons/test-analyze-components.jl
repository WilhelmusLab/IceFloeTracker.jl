using Images
using ImageComponentAnalysis
using DelimitedFiles
test_img = load("bw_new.tif")

# First label_components in test_img using 8-connectivity
@time Labels = label_components(test_img, trues(3,3),);

# Create measurements object
basic = BasicMeasurement(area = true, perimeter=false)
# elliptical = EllipseRegion(centroid=true,
#                            semiaxes=true,)

# Warning! took ~940 seconds to run on Carlos' system
# Output of line 26 follows: 
# 940.401786 seconds (85 allocations: 276.790 MiB, 0.00% gc time)
# 2479×8 DataFrame
#   Row │ l      Q₀        Q₁     Q₂     Q₃     Q₄       Qₓ     area
#       │ Int64  Int64     Int64  Int64  Int64  Int64    Int64  Float64
# ──────┼─────────────────────────────────────────────────────────────────────
#     1 │     1  36158671    797   2374    837    74016      8  76140.6
#     2 │     2  36232394     59    216     55     3979      0   4149.88
#   ⋮   │   ⋮       ⋮        ⋮      ⋮      ⋮       ⋮       ⋮          ⋮
#  2479 │  2479  36236651     12      8      8       24      0     38.0
#                                                            2476 rows omitted
@time components = analyze_components(Labels,basic)
writedlm("jl_imagecomp_areas.csv", components.area, ",")

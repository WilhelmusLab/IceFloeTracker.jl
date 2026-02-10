import ImageBinarization: AbstractImageBinarizationAlgorithm, binarize, AdaptiveThreshold
import Images:
    build_histogram,
    Colorant,
    AbstractGray,
    AbstractRGB,
    TransparentRGB,
    RGB,
    Gray,
    TransparentGray,
    red,
    green,
    blue,
    alpha,
    alphacolor
import Peaks: findmaxima, peakproms!, peakwidths!
import DataFrames: DataFrame, sort, Not
import ..ImageUtils: masker

"""
Given the edges and counts from build_histogram, identify local maxima and return the location of the
largest local maximum that is bright enough that it is possibly sea ice. Locations are determined by 
the edges, which by default are the left bin edges. Note also that peaks defaults to the left side of
plateaus. Returns Inf if there are no non-zero parts of the histogram with bins larger than the possible
ice threshold, or if there are no detected peaks larger than the minimum prominence.
"""
function get_ice_peaks(
    edges,
    counts;
    possible_ice_threshold::Float64=0.30,
    minimum_prominence::Float64=0.01,
    window_size::Int64=3,
)
    size(counts)
    counts = counts[1:end]
    normalizer = sum(counts[edges .> possible_ice_threshold])
    # Normalize the possible sea ice section of the histogram. 
    # Images with a lot of masked pixels can have large peaks at 0, which
    # we don't want to include in the normalization. If no potential
    # ice pixels, then return early
    counts = normalizer > 0 ? counts ./ normalizer : return Inf
    pks = findmaxima(counts, window_size) |> peakproms! |> peakwidths!
    pks_df = DataFrame(pks[Not(:data)])
    pks_df = sort(pks_df, :proms; rev=true)
    mx, argmx = findmax(pks_df.proms)
    mx < minimum_prominence && return Inf
    return edges[pks_df[argmx, :indices]]
end

"""
    IceDetectionAlgorithm

Functors to detect ice regions in an image.

Each algorithm `a` with parameters `kwargs...` can be called like:
- `binarize(image, a(; kwargs...))` 
- or `a(; kwargs...)(image)`.

"""
abstract type IceDetectionAlgorithm <: AbstractImageBinarizationAlgorithm end

function (a::IceDetectionAlgorithm)(image::AbstractArray{<:Colorant})
    return binarize(image, a)
end

"""
    IceDetectionThresholdMODIS721(;
        band_7_threshold::Real,
        band_2_threshold::Real,
        band_1_threshold::Real,
    )(image)
    binarize(
        modis_721_image, 
        a::IceDetectionThresholdMODIS721
    )
Returns pixels for a MODIS image where (band_7 < threshold AND band_2 > threshold AND band_1 > threshold).
"""
@kwdef struct IceDetectionThresholdMODIS721 <: IceDetectionAlgorithm
    band_7_max::Real
    band_2_min::Real
    band_1_min::Real
end

function (f::IceDetectionThresholdMODIS721)(out, modis_721_image, args...; kwargs...)
    band_7 = red.(modis_721_image)
    band_2 = green.(modis_721_image)
    band_1 = blue.(modis_721_image)
    mask_band_7 = band_7 .< f.band_7_max
    mask_band_2 = band_2 .> f.band_2_min
    mask_band_1 = band_1 .> f.band_1_min
    alpha_binary = alpha.(alphacolor.(modis_721_image)) .> 0.5
    @. out = mask_band_7 * mask_band_2 * mask_band_1 * alpha_binary
end

"""
    IceDetectionBrightnessPeaksMODIS721(;
        band_7_max::Real,
        possible_ice_threshold::Real
        nbins=64
        minimum_prominence=0.01
        window_size=3
    )(image)
    binarize(
        modis_721_image, 
        a::IceDetectionBrightnessPeaksMODIS721
    )

Uses the histogram of the band 1 and band 2 reflectance to determine thresholds for identifying
bright ice pixels (e.g., snow-covered floes or ice thicker than its surroundings). The algorithm 
builds a histogram using `nbins` bins, and finds the largest peak such that the peak brightness
is larger than the `possible_ice_threshold` and has prominence larger than `minimum_prominence`
using a comparison window of size `window_size` (see docs for `Peaks.findmaxima`). If `join_method` = "intersect", then select pixels where
the band 1 brightness is larger than the band 1 peak and the band 2 brightness is larger than the band
2 peak. Otherwise, if`join_method` = "union", select pixels where either band 1 or band 2 is brighter
than the threshold criteria. Finally, since clouds tend to have higher reflectance in band 7, mask
pixels with the band 7 brightness larger than `band_7_max`. It is designed to be used with MODIS
false color 7-2-1 imagery.
"""
@kwdef struct IceDetectionBrightnessPeaksMODIS721 <: IceDetectionAlgorithm
    band_7_max::Real
    possible_ice_threshold::Real
    nbins::Int64=64
    minimum_prominence::Float64=0.01
    window_size::Int64=3
    join_method="intersection"
end

function (f::IceDetectionBrightnessPeaksMODIS721)(out, modis_721_image, args...; kwargs...)
    band_7 = red.(modis_721_image)
    band_2 = green.(modis_721_image)
    band_1 = blue.(modis_721_image)

    alpha_binary = alpha.(alphacolor.(modis_721_image)) .> 0.5

    get_band_peak = function (band)
        return get_ice_peaks(
            build_histogram(band .* alpha_binary, f.nbins; minval=0, maxval=1)...;
            possible_ice_threshold=f.possible_ice_threshold,
            minimum_prominence=f.minimum_prominence,
            window_size=f.window_size
        )
    end

    band_2_peak = get_band_peak(band_2)
    band_1_peak = get_band_peak(band_1)

    mask_band_7 = band_7 .< f.band_7_max
    mask_band_2 = band_2 .> band_2_peak
    mask_band_1 = band_1 .> band_1_peak

    join_method = f.join_method
    if join_method ∉ ["intersection", "union"]
        @warn("Join method ", join_method, "not defined, defaulting to intersection")
        join_method = "intersection"
    end
    join_method == "intersection" && begin
        @. out = mask_band_7 && mask_band_2 && mask_band_1 && alpha_binary
    end
    join_method == "union" && begin
        @. out = mask_band_7 && (mask_band_2 || mask_band_1) && alpha_binary
    end
end

"""
    IceDetectionFirstNonZeroAlgorithm(;
        algorithms::Vector{IceDetectionAlgorithm},
    )(image)
    binarize(image, algorithms::IceDetectionFirstNonZeroAlgorithm)

Runs each algorithm from `algorithms` on the image, and returns the first which detects any ice.
"""
@kwdef struct IceDetectionFirstNonZeroAlgorithm <: IceDetectionAlgorithm
    algorithms::Vector{IceDetectionAlgorithm}
    threshold::Int64
end

function (f::IceDetectionFirstNonZeroAlgorithm)(out, img, args...; kwargs...)
    for algorithm in f.algorithms
        @debug algorithm
        result = binarize(img, algorithm)
        ice_sum = sum(result)
        if f.threshold < ice_sum
            @. out = result
            return nothing
        end
    end
    # In case we don't find anything, we're going to return zeros
    @. out = zero(eltype(out))
    return nothing
end

"""
    IceDetectionLopezAcosta2019(;
        band_7_threshold::Float64=Float64(5 / 255),
        band_2_threshold::Float64=Float64(230 / 255),
        band_1_threshold::Float64=Float64(240 / 255),
        band_7_threshold_relaxed::Float64=Float64(10 / 255),
        band_1_threshold_relaxed::Float64=Float64(190 / 255),
        possible_ice_threshold::Float64=Float64(75 / 255),
    )

Returns the first non-zero result of two threshold-based and one brightness-peak based ice detections.

Default thresholds are defined in the published Ice Floe Tracker article: Remote Sensing of the Environment 234 (2019) 111406.

"""
function IceDetectionLopezAcosta2019(;
    band_7_max::Float64=Float64(5 / 255),
    band_2_min::Float64=Float64(230 / 255),
    band_1_min::Float64=Float64(240 / 255),
    band_7_max_relaxed::Float64=Float64(10 / 255),
    band_1_min_relaxed::Float64=Float64(190 / 255),
    possible_ice_threshold::Float64=Float64(75 / 255),
)
    return IceDetectionFirstNonZeroAlgorithm([
        IceDetectionThresholdMODIS721(;
            band_7_max=band_7_max, band_2_min=band_2_min, band_1_min=band_1_min
        ),
        IceDetectionThresholdMODIS721(;
            band_7_max=band_7_max_relaxed,
            band_2_min=band_2_min,
            band_1_min=band_1_min_relaxed,
        ),
        IceDetectionBrightnessPeaksMODIS721(;
            band_7_max=band_7_max, possible_ice_threshold=possible_ice_threshold
        ),
    ], 1)
end

"""
    find_ice_labels(falsecolor_image, landmask; band_7_threshold, band_2_threshold, band_1_threshold, band_7_relaxed_threshold, band_1_relaxed_threshold, possible_ice_threshold)

Returns pixel indices of likely ice from false color reflectance image, using the thresholds from the Ice Floe Tracker article: Remote Sensing of the Environment 234 (2019) 111406.

# Arguments
- `falsecolor_image`: corrected reflectance false color image - bands [7,2,1]
- `landmask`: bitmatrix landmask for region of interest
- `band_7_threshold`: threshold value used to identify ice in band 7, N0f8(RGB intensity/255)
- `band_2_threshold`: threshold value used to identify ice in band 2, N0f8(RGB intensity/255)
- `band_1_threshold`: threshold value used to identify ice in band 2, N0f8(RGB intensity/255)
- `band_7_relaxed_threshold`: threshold value used to identify ice in band 7 if not found on first pass, N0f8(RGB intensity/255)
- `band_1_relaxed_threshold`: threshold value used to identify ice in band 1 if not found on first pass, N0f8(RGB intensity/255)

"""
function find_ice_labels(
    falsecolor_image::Matrix{RGB{Float64}}, landmask::BitMatrix; kwargs...
)
    ice_labels = find_ice_mask(falsecolor_image, landmask) |> get_ice_labels
    return ice_labels
end

function find_ice_mask(
    falsecolor_image::Matrix{RGB{Float64}}, landmask::BitMatrix; kwargs...
)
    masked_image = masker(landmask)(falsecolor_image)
    algorithm = IceDetectionLopezAcosta2019(; kwargs...)
    ice_mask = binarize(masked_image, algorithm)
    return ice_mask
end

function get_ice_labels(ice::AbstractArray{<:TransparentGray})
    return findall(vec(gray.(ice) .* alpha.(ice)) .> 0)
end

function get_ice_labels(ice::AbstractArray{<:AbstractGray})
    return findall(vec(gray.(ice)) .> 0)
end

# TODO: Find the right home for this function
# TODO: Set up wrapper for applying binarization to a tiled iterator. Let this
# method be one functor algorithm option, and the version that uses the Peaks as a second.
"""
tiled_adaptive_binarization(img, tiles; minimum_window_size=). 

    Applies the (AdaptiveThreshold)[https://juliaimages.org/ImageBinarization.jl/v0.1/#Adaptive-Threshold-1] binarization algorithm
    to each tile in the image. Following the recommendations from ImageBinarization, the default is to use the integer window size
    nearest to 1/8th the tile size if the tile is large enough. So that the window is large enough to include moderately large floes,
    the default minimum window size is 100 pixels (25 km for MODIS imagery). The minimum brightness parameter masks pixels with low
    grayscale intensity to prevent dark regions from getting brightened (i.e., the center of a large patch of open water).
    The "threshold_percentage" parameter is passed to the the AdaptiveThreshold function (percentage parameter).
"""
function tiled_adaptive_binarization(
    img, tiles; minimum_window_size=50, minimum_brightness=75 / 255, threshold_percentage=15
)
    canvas = zeros(size(img))
    img = deepcopy(img)
    img[Gray.(img) .< minimum_brightness] .= 0
    for tile in tiles
        L = Int64.(round.(minimum(length.(tile)) / 8, digits=0))
        L < minimum_window_size && (L = minimum_window_size)

        f = AdaptiveThreshold(img[tile...]; window_size=L, percentage=threshold_percentage)
        canvas[tile...] = binarize(img[tile...], f)
    end

    canvas[Gray.(img) .< minimum_brightness] .= 0
    return canvas
end

"""
    IceDetectionAlgorithm

Functors to detect ice regions in an image.

Each algorithm `a` with parameters `kwargs...` can be called like:
- `find_ice(image, a(; kwargs...))` 
- or `a(; kwargs...)(image)`.

"""
abstract type IceDetectionAlgorithm end

function (a::IceDetectionAlgorithm)(img; kwargs...)
    return find_ice(img, a; kwargs...)
end

"""
    IceDetectionThresholdMODIS721(;
        band_7_threshold::Real,
        band_2_threshold::Real,
        band_1_threshold::Real,
    )(image)
    find_ice(
        modis_721_image, 
        a::IceDetectionThresholdMODIS721
    )

Returns pixels for a MODIS image where (band_7 < threshold AND band_2 > threshold AND band_1 > threshold).
"""
@kwdef struct IceDetectionThresholdMODIS721 <: IceDetectionAlgorithm
    band_7_threshold::Real
    band_2_threshold::Real
    band_1_threshold::Real
end

function find_ice(
    modis_721_image::AbstractArray{<:Union{AbstractRGB,TransparentRGB}},
    a::IceDetectionThresholdMODIS721,
)
    band_7 = red.(modis_721_image)
    band_2 = green.(modis_721_image)
    band_1 = blue.(modis_721_image)

    mask_ice_band_7 = band_7 .< a.band_7_threshold
    mask_ice_band_2 = band_2 .> a.band_2_threshold
    mask_ice_band_1 = band_1 .> a.band_1_threshold

    ice = (mask_ice_band_7 .* mask_ice_band_2 .* mask_ice_band_1)

    ice_img = coloralpha.(Gray.(N0f8.(ice)), alpha.(alphacolor.(modis_721_image)))

    return ice_img
end

"""
    IceDetectionBrightnessPeaksMODIS721(;
        band_7_threshold::Real,
        possible_ice_threshold::Real
    )(image)
    find_ice(
        modis_721_image, 
        a::IceDetectionBrightnessPeaksMODIS721
    )

Returns pixels for a MODIS image where (band_7 < threshold AND both (band_2, band_1) are are above a peak value above some threshold).
"""
@kwdef struct IceDetectionBrightnessPeaksMODIS721 <: IceDetectionAlgorithm
    band_7_threshold::Real
    possible_ice_threshold::Real
end

function find_ice(
    modis_721_image::AbstractArray{<:Union{AbstractRGB,TransparentRGB}},
    a::IceDetectionBrightnessPeaksMODIS721,
)
    alpha_binary = alpha.(alphacolor.(modis_721_image)) .> 0.5
    band_7 = red.(modis_721_image)
    band_2 = green.(modis_721_image)
    band_1 = blue.(modis_721_image)

    mask_ice_band_7 = band_7 .< a.band_7_threshold
    band_2_peak = _find_reflectance_peaks(band_2 .* alpha_binary; a.possible_ice_threshold)
    band_1_peak = _find_reflectance_peaks(band_1 .* alpha_binary; a.possible_ice_threshold)

    mask_ice_band_2 = band_2 .> band_2_peak
    mask_ice_band_1 = band_1 .> band_1_peak
    ice = mask_ice_band_7 .* mask_ice_band_2 .* mask_ice_band_1 .* alpha_binary

    ice_img = coloralpha.(Gray.(N0f8.(ice)), alpha.(alphacolor.(modis_721_image)))

    return ice_img
end

function _find_reflectance_peaks(
    reflectance_channel::AbstractArray{<:Real}; possible_ice_threshold::Real=N0f8(75 / 255)
)
    reflectance_channel[reflectance_channel .< possible_ice_threshold] .= 0
    _, counts = ImageContrastAdjustment.build_histogram(reflectance_channel)
    locs, _ = Peaks.findmaxima(counts)
    sort!(locs; rev=true)
    return locs[2] / 255.0 # second greatest peak
end

"""
    IceDetectionFirstNonZeroAlgorithm(;
        algorithms::Vector{IceDetectionAlgorithm},
    )(image)
    find_ice(image, algorithms::IceDetectionFirstNonZeroAlgorithm)

Runs each algorithm from `algorithms` on the image, and returns the first which detects any ice.
"""
@kwdef struct IceDetectionFirstNonZeroAlgorithm <: IceDetectionAlgorithm
    algorithms::Vector{IceDetectionAlgorithm}
end

function find_ice(
    image::AbstractArray{<:Union{AbstractRGB,TransparentRGB}},
    a::IceDetectionFirstNonZeroAlgorithm,
)
    let ice
        for algorithm in a.algorithms
            @debug algorithm
            ice = find_ice(image, algorithm)
            ice_sum = sum(gray.(ice) .* alpha.(ice))
            if ice_sum > 0
                break
            end
        end
        return ice
    end
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
    band_7_threshold::Float64=Float64(5 / 255),
    band_2_threshold::Float64=Float64(230 / 255),
    band_1_threshold::Float64=Float64(240 / 255),
    band_7_threshold_relaxed::Float64=Float64(10 / 255),
    band_1_threshold_relaxed::Float64=Float64(190 / 255),
    possible_ice_threshold::Float64=Float64(75 / 255),
)
    return IceDetectionFirstNonZeroAlgorithm([
        IceDetectionThresholdMODIS721(;
            band_7_threshold=band_7_threshold,
            band_2_threshold=band_2_threshold,
            band_1_threshold=band_1_threshold,
        ),
        IceDetectionThresholdMODIS721(;
            band_7_threshold=band_7_threshold_relaxed,
            band_2_threshold=band_2_threshold,
            band_1_threshold=band_1_threshold_relaxed,
        ),
        IceDetectionBrightnessPeaksMODIS721(;
            band_7_threshold=band_7_threshold, possible_ice_threshold=possible_ice_threshold
        ),
    ])
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
)::Vector{Int64}
    masked_image = masker(.!(landmask))(falsecolor_image)
    algorithm = IceDetectionLopezAcosta2019(; kwargs...)
    ice = IceFloeTracker.find_ice(masked_image, algorithm)
    ice_labels = get_ice_labels(ice)
    return ice_labels
end

function get_ice_labels(ice::AbstractArray{<:TransparentGray})
    return findall(vec(gray.(ice) .* alpha.(ice)) .> 0)
end

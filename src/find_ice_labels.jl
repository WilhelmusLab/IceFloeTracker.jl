# Type of ref_image_band_2
__T__ = SubArray{
    Float64,
    2,
    Base.ReinterpretArray{Float64,3,RGB{Float64},Matrix{RGB{Float64}},true},
    Tuple{Int64,Base.Slice{Base.OneTo{Int64}},Base.Slice{Base.OneTo{Int64}}},
    false,
}

"""
    find_reflectance_peaks(reflectance_channel, possible_ice_threshold;)

Find histogram peaks in single channels of a reflectance image and return the second greatest peak. If needed, edges can be returned as the first object from `build_histogram`. Similarly, peak values can be returned as the second object from `findmaxima`.

# Arguments
- `reflectance_channel`: either band 2 or band 1 of false-color reflectance image
- `possible_ice_threshold`: threshold value used to identify ice if not found on first or second pass

"""
function find_reflectance_peaks(
    reflectance_channel::Union{__T__,Matrix{Float64}};
    possible_ice_threshold::Float64=Float64(75 / 255),
)::Int64
    reflectance_channel[reflectance_channel .< possible_ice_threshold] .= 0 #75 / 255
    _, counts = ImageContrastAdjustment.build_histogram(reflectance_channel)
    locs, _ = Peaks.findmaxima(counts)
    sort!(locs; rev=true)
    return locs[2] # second greatest peak
end

"""
    find_ice_labels(falsecolor_image, landmask; band_7_threshold, band_2_threshold, band_1_threshold, band_7_relaxed_threshold, band_1_relaxed_threshold, possible_ice_threshold)

Locate the pixels of likely ice from false color reflectance image. Returns a binary mask with ice floes contrasted from background. Default thresholds are defined in the published Ice Floe Tracker article: Remote Sensing of the Environment 234 (2019) 111406.

# Arguments
- `modis_721`: corrected reflectance false color image - bands [7,2,1]
- `landmask`: bitmatrix landmask for region of interest
- `band_7_threshold`: threshold value used to identify ice in band 7, N0f8(RGB intensity/255)
- `band_2_threshold`: threshold value used to identify ice in band 2, N0f8(RGB intensity/255)
- `band_1_threshold`: threshold value used to identify ice in band 2, N0f8(RGB intensity/255)
- `band_7_relaxed_threshold`: threshold value used to identify ice in band 7 if not found on first pass, N0f8(RGB intensity/255)
- `band_1_relaxed_threshold`: threshold value used to identify ice in band 1 if not found on first pass, N0f8(RGB intensity/255)

"""
function find_ice_labels(
    modis_721::AbstractArray{T},
    not_land::BitMatrix;
    band_7_threshold::Real=(5 / 255),
    band_2_threshold::Real=(230 / 255),
    band_1_threshold::Real=(240 / 255),
    band_7_threshold_relaxed::Real=(10 / 255),
    band_1_threshold_relaxed::Real=(190 / 255),
    possible_ice_threshold::Real=(75 / 255),
)::Vector{Int64} where {T<:Union{AbstractRGB,TransparentRGB}}
    modis_band07 = red.(modis_721)
    modis_band02 = green.(modis_721)
    modis_band01 = blue.(modis_721)

    mask_ice_band_7 = modis_band07 .< band_7_threshold
    mask_ice_band_2 = modis_band02 .> band_2_threshold
    mask_ice_band_1 = modis_band01 .> band_1_threshold
    ice = mask_ice_band_7 .* mask_ice_band_2 .* mask_ice_band_1
    ice_labels = remove_landmask(not_land, ice)

    if sum(abs.(ice_labels)) != 0
        return ice_labels
    end

    mask_ice_band_7 = modis_band07 .< band_7_threshold_relaxed #10 / 255
    mask_ice_band_1 = modis_band02 .> band_1_threshold_relaxed #190 / 255
    ice = mask_ice_band_7 .* mask_ice_band_2 .* mask_ice_band_1
    ice_labels = remove_landmask(not_land, ice)

    if sum(abs.(ice_labels)) != 0
        return ice_labels
    end

    band_2_peak = find_reflectance_peaks(
        modis_band02; possible_ice_threshold=possible_ice_threshold
    )
    band_1_peak = find_reflectance_peaks(
        modis_band01; possible_ice_threshold=possible_ice_threshold
    )
    mask_ice_band_2 = modis_band02 .> band_2_peak / 255
    mask_ice_band_1 = modis_band01 .> band_1_peak / 255
    ice = mask_ice_band_7 .* mask_ice_band_2 .* mask_ice_band_1
    ice_labels = remove_landmask(not_land, ice)

    return ice_labels
end

"""
    find_ice_labels(falsecolor_image, landmask; band_7_threshold, band_2_threshold, band_1_threshold, band_7_relaxed_threshold, band_1_relaxed_threshold, possible_ice_threshold)

Locate the pixels of likely ice from false color reflectance image. Returns a binary mask with ice floes contrasted from background. Default thresholds are defined in the published Ice Floe Tracker article: Remote Sensing of the Environment 234 (2019) 111406.

# Arguments
- `modis_721`: corrected reflectance false color image - bands [7,2,1]
- `landmask`: bitmatrix landmask for region of interest
- `band_7_threshold`: threshold value used to identify ice in band 7, N0f8(RGB intensity/255)
- `band_2_threshold`: threshold value used to identify ice in band 2, N0f8(RGB intensity/255)
- `band_1_threshold`: threshold value used to identify ice in band 2, N0f8(RGB intensity/255)
- `band_7_relaxed_threshold`: threshold value used to identify ice in band 7 if not found on first pass, N0f8(RGB intensity/255)
- `band_1_relaxed_threshold`: threshold value used to identify ice in band 1 if not found on first pass, N0f8(RGB intensity/255)

"""
function find_ice(
    modis_721::AbstractArray{T},
    not_land::AbstractArray{<:Gray};
    band_7_threshold::Real=(5 / 255),
    band_2_threshold::Real=(230 / 255),
    band_1_threshold::Real=(240 / 255),
    band_7_threshold_relaxed::Real=(10 / 255),
    band_1_threshold_relaxed::Real=(190 / 255),
    possible_ice_threshold::Real=(75 / 255),
) where {T<:Union{AbstractRGB,TransparentRGB}}
    modis_band07 = red.(modis_721)
    modis_band02 = green.(modis_721)
    modis_band01 = blue.(modis_721)

    mask_ice_band_7 = modis_band07 .< band_7_threshold
    mask_ice_band_2 = modis_band02 .> band_2_threshold
    mask_ice_band_1 = modis_band01 .> band_1_threshold
    ice = mask_ice_band_7 .* mask_ice_band_2 .* mask_ice_band_1
    ice = IceFloeTracker.apply_landmask(ice, not_land)

    if count(ice) != 0  # ice is a gray image – but we need something which is an image but also a boolean
        return Gray.(ice)
    end

    mask_ice_band_7 = modis_band07 .< band_7_threshold_relaxed
    mask_ice_band_2 = modis_band02 .> band_2_threshold
    mask_ice_band_1 = modis_band02 .> band_1_threshold_relaxed
    ice = mask_ice_band_7 .* mask_ice_band_2 .* mask_ice_band_1
    ice = IceFloeTracker.apply_landmask(ice, not_land)

    if count(ice) != 0
        return Gray.(ice)
    end

    mask_ice_band_7 = modis_band07 .< band_7_threshold_relaxed
    band_2_peak = find_reflectance_peaks(
        modis_band02; possible_ice_threshold=possible_ice_threshold
    )
    mask_ice_band_2 = modis_band02 .> band_2_peak / 255
    band_1_peak = find_reflectance_peaks(
        modis_band01; possible_ice_threshold=possible_ice_threshold
    )
    mask_ice_band_1 = modis_band01 .> band_1_peak / 255
    ice = mask_ice_band_7 .* mask_ice_band_2 .* mask_ice_band_1
    ice = IceFloeTracker.apply_landmask(ice, not_land)

    return Gray.(ice)
end

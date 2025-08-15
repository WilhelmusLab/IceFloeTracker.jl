
function find_reflectance_peaks(
    reflectance_channel::AbstractArray{<:Real}; possible_ice_threshold::Real=N0f8(75 / 255)
)
    reflectance_channel[reflectance_channel .< possible_ice_threshold] .= 0
    _, counts = ImageContrastAdjustment.build_histogram(reflectance_channel)
    locs, _ = Peaks.findmaxima(counts)
    sort!(locs; rev=true)
    return locs[2] / 255.0 # second greatest peak
end

function get_ice_labels(ice::AbstractArray{<:TransparentGray})
    return findall(vec(gray.(ice) .* alpha.(ice)) .> 0)
end

"""
    find_ice_labels(falsecolor_image, landmask; band_7_threshold, band_2_threshold, band_1_threshold, band_7_relaxed_threshold, band_1_relaxed_threshold, possible_ice_threshold)

Locate the pixels of likely ice from false color reflectance image. Returns a binary mask with ice floes contrasted from background. Default thresholds are defined in the published Ice Floe Tracker article: Remote Sensing of the Environment 234 (2019) 111406.

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
    falsecolor_image::Matrix{RGB{Float64}},
    landmask::BitMatrix;
    band_7_threshold::Float64=Float64(5 / 255),
    band_2_threshold::Float64=Float64(230 / 255),
    band_1_threshold::Float64=Float64(240 / 255),
    band_7_threshold_relaxed::Float64=Float64(10 / 255),
    band_1_threshold_relaxed::Float64=Float64(190 / 255),
    possible_ice_threshold::Float64=Float64(75 / 255),
)::Vector{Int64}

    ## Make ice masks
    cv = channelview(falsecolor_image)

    @info "first attempt at finding ice labels"
    mask_ice_band_7 = @view(cv[1, :, :]) .< band_7_threshold #5 / 255
    mask_ice_band_2 = @view(cv[2, :, :]) .> band_2_threshold #230 / 255
    mask_ice_band_1 = @view(cv[3, :, :]) .> band_1_threshold #240 / 255
    ice = mask_ice_band_7 .* mask_ice_band_2 .* mask_ice_band_1
    ice_labels = apply_landmask(ice, landmask; as_indices=true)
    # @info "Done with masks" # to uncomment when logger is added

    ## Find likely ice floes
    if sum(abs.(ice_labels)) == 0
        @info "second attempt at finding ice labels"
        mask_ice_band_7 = @view(cv[1, :, :]) .< band_7_threshold_relaxed #10 / 255
        mask_ice_band_1 = @view(cv[3, :, :]) .> band_1_threshold_relaxed #190 / 255
        ice = mask_ice_band_7 .* mask_ice_band_2 .* mask_ice_band_1
        ice_labels = apply_landmask(ice, landmask; as_indices=true)
        if sum(abs.(ice_labels)) == 0
            @info "third attempt at finding ice labels"
            ref_image_band_2 = @view(cv[2, :, :])
            ref_image_band_1 = @view(cv[3, :, :])
            band_2_peak = find_reflectance_peaks(
                ref_image_band_2; possible_ice_threshold=possible_ice_threshold
            )
            band_1_peak = find_reflectance_peaks(
                ref_image_band_1; possible_ice_threshold=possible_ice_threshold
            )
            mask_ice_band_2 = @view(cv[2, :, :]) .> band_2_peak / 255
            mask_ice_band_1 = @view(cv[3, :, :]) .> band_1_peak / 255
            ice = mask_ice_band_7 .* mask_ice_band_2 .* mask_ice_band_1
            ice_labels = apply_landmask(ice, landmask; as_indices=true)
        end
    end
    # @info "Done with ice labels" # to uncomment when logger is added
    return ice_labels
end

abstract type IceDetectionAlgorithm end

function (a::IceDetectionAlgorithm)(img; kwargs...)
    return find_ice(img, a; kwargs...)
end

@kwdef struct IceDetectionFirstNonZeroAlgorithm <: IceDetectionAlgorithm
    algorithms::Vector{IceDetectionAlgorithm}
end

function find_ice(
    modis_721_image::AbstractArray{<:Union{AbstractRGB,TransparentRGB}},
    a::IceDetectionFirstNonZeroAlgorithm,
)
    let ice
        for algorithm in a.algorithms
            @debug algorithm
            ice = find_ice(modis_721_image, algorithm)
            ice_sum = sum(gray.(ice) .* alpha.(ice))
            if ice_sum > 0
                break
            end
        end
        return ice
    end
end

@kwdef struct IceDetectionThreshold <: IceDetectionAlgorithm
    band_7_threshold::Real
    band_2_threshold::Real
    band_1_threshold::Real
end

function find_ice(
    modis_721_image::AbstractArray{<:Union{AbstractRGB,TransparentRGB}},
    a::IceDetectionThreshold,
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

@kwdef struct IceDetectionBrightnessPeaks <: IceDetectionAlgorithm
    band_7_threshold::Real
    possible_ice_threshold::Real
end

function find_ice(
    modis_721_image::AbstractArray{<:Union{AbstractRGB,TransparentRGB}},
    a::IceDetectionBrightnessPeaks,
)
    alpha_binary = alpha.(alphacolor.(modis_721_image)) .> 0.5
    band_7 = red.(modis_721_image)
    band_2 = green.(modis_721_image)
    band_1 = blue.(modis_721_image)

    mask_ice_band_7 = band_7 .< a.band_7_threshold
    band_2_peak = find_reflectance_peaks(band_2 .* alpha_binary; a.possible_ice_threshold)
    band_1_peak = find_reflectance_peaks(band_1 .* alpha_binary; a.possible_ice_threshold)

    mask_ice_band_2 = band_2 .> band_2_peak
    mask_ice_band_1 = band_1 .> band_1_peak
    ice = mask_ice_band_7 .* mask_ice_band_2 .* mask_ice_band_1 .* alpha_binary

    ice_img = coloralpha.(Gray.(N0f8.(ice)), alpha.(alphacolor.(modis_721_image)))

    return ice_img
end

function LopezAcosta2019IceDetection()
    return IceDetectionFirstNonZeroAlgorithm([
        IceDetectionThreshold(;
            band_7_threshold=N0f8(5 / 255),
            band_2_threshold=N0f8(230 / 255),
            band_1_threshold=N0f8(240 / 255),
        ),
        IceDetectionThreshold(;
            band_7_threshold=N0f8(10 / 255),
            band_2_threshold=N0f8(230 / 255),
            band_1_threshold=N0f8(190 / 255),
        ),
        IceDetectionBrightnessPeaks(;
            band_7_threshold=N0f8(5 / 255), possible_ice_threshold=N0f8(75 / 255)
        ),
    ])
end
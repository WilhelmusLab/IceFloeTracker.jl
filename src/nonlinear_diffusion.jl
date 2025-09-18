"""Nonlinear Diffusion

Anisotropic diffusion was introduced by Perona and Malik (1987) and refined in subsequent publications.
We include two implementations of the Perona and Malik algorithm: one adapted from https://github.com/Red-Portal/UltrasoundDesignGallery.jl, and the
other a reimplementation of approach used in the default Matlab image processing toolbox.
Since the algorithm is not truly anisotropic, we refer to it instead as nonlinear diffusion.

PeronaMalikDiffusion(λ, K, niters, g)
    λ = Parameter weighting the diffusion rate, needs to be between 0 and 0.25 for stability.
    K = Numerator for the image gradient function. (TBD: Option to estimate from image gradient histogram)
    niters = Number of interations
    g = "exponential", "inverse_quadratic" (TBD: Option to provide user-defined function)


P. Perona and J. Malik (November 1987). "Scale-space and edge detection using anisotropic diffusion". Proceedings of IEEE Computer Society Workshop on Computer Vision. pp. 16–22.
P. Perona and J. Malik, Scale-Space and Edge Detection Using Anisotropic Diffusion, IEEE Transactions on Pattern Analysis and Machine Intelligence, 12(7):629-639, July 1990
G. Grieg, O. Kubler, R. Kikinis, and F. A. Jolesz, Nonlinear Anisotropic Filtering of MRI Data, IEEE Transactions on Medical Imaging, 11(2):221-232, June 1992

"""

include("gradient_functions.jl")
abstract type AbstractDiffusionAlgorithm end

@kwdef struct PeronaMalikDiffusion <: AbstractDiffusionAlgorithm
    λ::Float64 = 0.1
    K::Number = 0.1
    niters::Int = 3
    g::String = "inverse_quadratic"

    # enforce conditions
    function PeronaMalikDiffusion(λ, K, niters, g)
        !(0 < λ <= 0.25) && throw(ArgumentError("Lambda must be between zero and 0.25"))
        K <= 0 && throw(ArgumentError("K must be greater than zero"))
        niters <= 0 &&
            throw(ArgumentError("Number of iterations must be greater than zero"))
        g ∉ SUPPORTED_GRADIENT_FUNCTIONS && throw(
            ArgumentError(
                "Unknown function name. Supported functions: $SUPPORTED_GRADIENT_FUNCTIONS",
            ),
        )
        return new(λ, K, niters, g)
    end
end

# Default to using Perona Malik diffusion. Future releases may include more modern algorithms.
function nonlinear_diffusion(
    img::AbstractArray{<:Union{AbstractRGB,TransparentRGB,AbstractGray}},
    f::AbstractDiffusionAlgorithm=PeronaMalikDiffusion(),
)
    return f(img)
end

function nonlinear_diffusion(
    img::AbstractArray{<:Union{AbstractRGB,TransparentRGB,AbstractGray}},
    λ::Float64,
    K::Number,
    niters::Int,
)
    return nonlinear_diffusion(img, PeronaMalikDiffusion(λ, K, niters, "inverse_quadratic"))
end

function (f::PeronaMalikDiffusion)(img::AbstractArray{<:AbstractGray})
    # Get the gradient function from the supported functions
    g = SUPPORTED_GRADIENT_FUNCTIONS[f.g]

    # Future option: Implement updater using 8-connectivity instead of 4-connectivity
    function pmd_updater!(image, output, g, λ, k)
        M, N = size(image)
        padded_array = padarray(image, Pad(:replicate, 1, 1))

        ∇n = padded_array[1:M, 0:(N - 1)] .- image
        ∇s = padded_array[1:M, 2:(N + 1)] .- image
        ∇e = padded_array[2:(M + 1), 1:N] .- image
        ∇w = padded_array[0:(M - 1), 1:N] .- image

        Cn = g.(abs.(∇n), k)
        Cs = g.(abs.(∇s), k)
        Ce = g.(abs.(∇e), k)
        Cw = g.(abs.(∇w), k)

        return output .= image + λ .* (Cn .* ∇n + Cs .* ∇s + Ce .* ∇e + Cw .* ∇w)
    end

    # Since we are doing math on the image, we need to reinterpret as float
    _img = float64.(img)
    _out = deepcopy(_img)

    for _ in 1:(f.niters)
        # Future option: estimate k from the data, as in P-M paper
        pmd_updater!(_img, _out, g, f.λ, f.K)
        _img, _out = _out, _img
    end

    # Map the diffused image back to the original base color type
    recast_img_type = base_color_type(eltype(_out)){eltype(eltype(img))}
    return recast_img_type.(_out)
end

function (f::PeronaMalikDiffusion)(img::AbstractArray{<:Union{AbstractRGB,TransparentRGB}})
    # TBD: loop through colorview applying the diffusion function
    cv = channelview(img)
    for i in 1:3
        cvi_gray = Gray.(cv[i, :, :])
        diffused_cvi = f(cvi_gray)
        cv[i, :, :] .= Float64.(diffused_cvi)
    end

    return colorview(eltype(img), cv)
    # TBD
end

# We can replace this with a loop through the channel view, as done for the function above
function anisotropic_diffusion_3D(I)
    rgbchannels = get_rgb_channels(I)

    for i in 1:3
        rgbchannels[:, :, i] .= anisotropic_diffusion_2D(rgbchannels[:, :, i])
    end

    return rgbchannels
end

function anisotropic_diffusion_2D(
    # Implementation of the matlab 2D anisotropic diffusion filter default mode
    # by Carlos Paniagua
    I::AbstractMatrix{T};
    gradient_threshold::Union{T,Nothing}=nothing,
    niter::Int=1,
) where {T}
    if eltype(I) <: Int
        I = Gray.(I ./ 255)
    end

    # Determine the gradient threshold if not provided
    # dmw: this is more of a scaling factor than a threshold is it not?
    if gradient_threshold === nothing
        dynamic_range = maximum(I) - minimum(I)
        gradient_threshold = 0.1 * dynamic_range
    end

    # Padding the image (corrected)
    padded_img = padarray(I, Pad(:replicate, (1, 1)))
    dd = sqrt(2)
    diffusion_rate = 1 / 8  # Fixed for maximal connectivity (8 neighbors)

    for _ in 1:niter
        # These are zero-indexed offset arrays
        diff_img_north =
            padded_img[0:(end - 1), 1:(end - 1)] .- padded_img[1:end, 1:(end - 1)]
        diff_img_east =
            padded_img[1:(end - 1), 1:end] .- padded_img[1:(end - 1), 0:(end - 1)]
        diff_img_nw = padded_img[0:(end - 2), 0:(end - 2)] .- I
        diff_img_ne = padded_img[0:(end - 2), 2:end] .- I
        diff_img_sw = padded_img[2:end, 0:(end - 2)] .- I
        diff_img_se = padded_img[2:end, 2:end] .- I

        # Exponential conduction coefficients
        conduct_coeff_north = exp.(-(abs.(diff_img_north) ./ gradient_threshold) .^ 2)
        conduct_coeff_east = exp.(-(abs.(diff_img_east) ./ gradient_threshold) .^ 2)
        conduct_coeff_nw = exp.(-(abs.(diff_img_nw) ./ gradient_threshold) .^ 2)
        conduct_coeff_ne = exp.(-(abs.(diff_img_ne) ./ gradient_threshold) .^ 2)
        conduct_coeff_sw = exp.(-(abs.(diff_img_sw) ./ gradient_threshold) .^ 2)
        conduct_coeff_se = exp.(-(abs.(diff_img_se) ./ gradient_threshold) .^ 2)

        # Flux calculations
        flux_north = conduct_coeff_north .* diff_img_north
        flux_east = conduct_coeff_east .* diff_img_east
        flux_nw = conduct_coeff_nw .* diff_img_nw
        flux_ne = conduct_coeff_ne .* diff_img_ne
        flux_sw = conduct_coeff_sw .* diff_img_sw
        flux_se = conduct_coeff_se .* diff_img_se

        # Back to regular 1-indexed arrays
        flux_north_diff = flux_north[1:(end - 1), :] .- flux_north[2:end, :]
        flux_east_diff = flux_east[:, 2:end] .- flux_east[:, 1:(end - 1)]

        # Discrete PDE solution
        sum_ = @. (1 / (dd^2)) * (flux_nw + flux_ne + flux_sw + flux_se)

        I = @. I + diffusion_rate * (flux_north_diff + flux_east_diff + sum_)
    end

    return I
end

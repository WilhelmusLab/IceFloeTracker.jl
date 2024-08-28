# Anisotropic Image Diffusion
# Adapted from https://github.com/Red-Portal/UltrasoundDesignGallery.jl
## MIT license with permission to use 

function pmad_kernel!(image, output, g, λ)
    M, N = size(image)

    @inbounds for j in 1:N
        @simd for i in 1:M
            w = image[max(i - 1, 1), j]
            n = image[i, max(j - 1, 1)]
            c = image[i, j]
            s = image[i, min(j + 1, N)]
            e = image[min(i + 1, M), j]

            ∇n = n - c
            ∇s = s - c
            ∇w = w - c
            ∇e = e - c

            Cn = g(abs(∇n))
            Cs = g(abs(∇s))
            Cw = g(abs(∇w))
            Ce = g(abs(∇e))

            output[i, j] = c + λ .* (Cn ⊙ ∇n + Cs ⊙ ∇s + Ce ⊙ ∇e + Cw ⊙ ∇w)
        end
    end
end

function invert_color(color::RGB{T}) where {T<:AbstractFloat}
    return RGB(1.0 / color.r, 1.0 / color.g, 1.0 / color.b)
end

function invert_color(color::Gray{T}) where {T<:AbstractFloat}
    return Gray(1.0 / color.val)
end

function diffusion(
    image::Matrix{T}, λ::Float64, K::Int, niters::Int
) where {T<:Color{Float64}}
    #=
        Perona, Pietro, and Jitendra Malik. 
        "Scale-space and edge detection using anisotropic diffusion." 
        IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI), 1990.
    =#
    if !(0 <= λ <= 0.25)
        error("Lambda must be between zero and 0.25")
    end

    @inline function g(norm∇I)
        coef = norm∇I / K
        denom = T(1) .+ coef ⊙ coef
        return invert_color(denom)
    end

    output = deepcopy(image)
    for _ in 1:niters
        pmad_kernel!(image, output, g, λ)
        image, output = output, image
    end
    return output
end

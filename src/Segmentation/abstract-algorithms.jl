using Images: AbstractRGB, TransparentRGB, Gray, float64
using ..Preprocessing: create_landmask

abstract type IceFloeSegmentationAlgorithm end

function (p::IceFloeSegmentationAlgorithm)(
    truecolor::T₁,
    falsecolor::T₂,
    landmask::T₃;
    intermediate_results_callback::Union{Nothing,Function}=nothing,
) where {
    T₁<:AbstractMatrix{<:Union{AbstractRGB,TransparentRGB}},
    T₂<:AbstractMatrix{<:Union{AbstractRGB,TransparentRGB}},
    T₃<:AbstractMatrix{<:Union{Bool,Gray,AbstractRGB,TransparentRGB}},
}
    landmask, coastal_buffer_mask = create_landmask(
        float64.(landmask), p.coastal_buffer_structuring_element
    )

    landmask = reinterpret(Bool, landmask)
    coastal_buffer_mask = reinterpret(Bool, coastal_buffer_mask)

    return p(
        truecolor,
        falsecolor,
        landmask,
        coastal_buffer_mask;
        intermediate_results_callback=intermediate_results_callback,
    )
end

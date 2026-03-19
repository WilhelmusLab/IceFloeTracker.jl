abstract type IceFloeSegmentationAlgorithm end

using Images: AbstractRGB, TransparentRGB, Gray, float64
using ..Preprocessing: create_landmask

function (p::IceFloeSegmentationAlgorithm)(
    truecolor::T₁,
    falsecolor::T₂,
    landmask::T₃;
    intermediate_results_callback::Union{Nothing,Function}=nothing,
) where {
    T₁<:AbstractMatrix{<:Union{AbstractRGB,TransparentRGB}},
    T₂<:AbstractMatrix{<:Union{AbstractRGB,TransparentRGB}},
    T₃<:AbstractMatrix{<:Union{Bool,Gray{Bool},AbstractRGB,TransparentRGB}},
}
    landmask, coastal_buffer_mask = create_landmask(
        float64.(landmask), p.coastal_buffer_structuring_element
    )
    return p(
        truecolor,
        falsecolor,
        landmask,
        coastal_buffer_mask;
        intermediate_results_callback=intermediate_results_callback,
    )
end
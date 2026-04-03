using Images: AbstractRGB, TransparentRGB, Gray, float64

abstract type IceFloeSegmentationAlgorithm end

import ..Preprocessing: create_coastal_buffer_mask

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
    landmask_ = reinterpret(Bool, landmask)
    coastal_buffer_mask = create_coastal_buffer_mask(
        landmask_, centered(p.coastal_buffer_structuring_element)
    )

    return p(
        truecolor,
        falsecolor,
        landmask,
        coastal_buffer_mask;
        intermediate_results_callback=intermediate_results_callback,
    )
end

import Images: SegmentedImage
import Dates: DateTime
import ..Segmentation: regionprops_table

abstract type AbstractTracker end

@kwdef struct FloeTracker <: AbstractTracker
    filter_function::AbstractFloeFilterFunction
    matching_function::AbstractFloeMatchingFunction
    # TODO: make area and time filters just part of the filter and matcher functions
    minimum_area::Real = 100
    maximum_area::Real = 90e3
    maximum_time_step::Period = Day(2)
end

function (t::FloeTracker)(
    segmented_images::Vector{<:SegmentedImage}, passtimes::Vector{DateTime}
)
    props = regionprops_table.(segmented_images)
    add_uuids!.(props)
    !issorted(passtimes) && @warn "Passtimes are not in ascending order."
    add_passtimes!.(props, passtimes)
    add_floemasks!.(props, segmented_images)
    add_ψs!.(props)

    tracking_results = floe_tracker(
        props,
        t.filter_function,
        t.matching_function;
        minimum_area=t.minimum_area,
        maximum_area=t.maximum_area,
        maximum_time_step=t.maximum_time_step,
    )

    return tracking_results
end
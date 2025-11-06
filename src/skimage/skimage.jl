module skimage

using PythonCall

export sk_morphology, sk_exposure, measure
const sk_morphology = PythonCall.pynew()
const sk_exposure = PythonCall.pynew()

function __init__()
    PythonCall.pycopy!(sk_morphology, pyimport("skimage.morphology"))
    PythonCall.pycopy!(sk_exposure, pyimport("skimage.exposure"))
    return nothing
end

module measure

    export regionprops_table

    using PythonCall
    import DataFrames: DataFrame

    const _measure = PythonCall.pynew()
    function __init__()
        PythonCall.pycopy!(_measure, PythonCall.pyimport("skimage.measure"))
        return nothing
    end

    """Wrapper around the python `skimage.measure.regionprops_table` function."""
    function regionprops_table(
        label_img::Matrix{<:Integer},
        intensity_img::Union{Nothing,AbstractMatrix}=nothing;
        properties::Union{Vector{<:AbstractString},Tuple{String,Vararg{String}}}=(
            "centroid",
            "area",
            "major_axis_length",
            "minor_axis_length",
            "convex_area",
            "bbox",
            "perimeter",
            "orientation",
        ),
    )::DataFrame
        result =
            _measure.regionprops_table(
                PythonCall.Py(label_img).to_numpy(),
                intensity_img ? PythonCall.Py(intensity_img).to_numpy() : nothing;
                properties,
            ) |>
            (x -> PythonCall.pyconvert(Dict, x)) |>
            DataFrame

        return result
    end
end

module morphology end

module exposure end

end

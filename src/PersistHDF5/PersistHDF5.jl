
module PersistHDF5

using HDF5, Images, Dates, TimeZones, DataFrames
import ..Geospatial: latlon
import ..Segmentation: regionprops_table, converttounits!

export make_hdf5

function choose_dtype(mx::T) where {T<:Integer}
    types = [UInt8, UInt16, UInt32, UInt64, UInt128]
    for t_ in types
        if typemin(t_) <= mx <= typemax(t_)
            return t_
        end
    end
    return error("$mx cannot be represented by any of $types")
end

function make_hdf5(
    output_path::AbstractString,
    ;
    passtime::ZonedDateTime,
    crs_ref_image_path::AbstractString,
    truecolor_path::AbstractString,
    falsecolor_path::AbstractString,
    labeled::AbstractMatrix,
    props::DataFrame,
    cloud_mask::AbstractMatrix,
    ice_mask::AbstractMatrix,
    landmask::AbstractMatrix,
    coastal_buffer_mask::AbstractMatrix,
    iftversion::VersionNumber=pkgversion(@__MODULE__),
    reference::AbstractString="https://doi.org/10.1016/j.rse.2019.111406",
    contact::AbstractString="mmwilhelmus@brown.edu",
)
    ptsunix = Int64(Dates.datetime2unix(DateTime(passtime)))
    latlondata = latlon(crs_ref_image_path)

    local crs_name::String
    try
        crs_name = Dict(
            3413 => "EPSG:3413 NSIDC north polar stereographic",
            3031 => "EPSG:3031 NSIDC south polar stereographic",
            4326 => "EPSG:4326 WGS84 lat/lon",
            3857 => "EPSG:3857 Web Mercator",
        )[latlondata[:crs]]
    catch e
        if e isa KeyError
            @warn "CRS $(latlondata[:crs]) not recognized. CRS will be recorded as EPSG:$(string(latlondata[:crs])) in the output file attributes, but no short name will be provided."
            crs_name = "EPSG:$(string(latlondata[:crs]))"
        else
            rethrow(e)
        end
    end

    h5open(output_path, "w") do file
        @info "Add top-level attributes"
        attrs(file)["fname_falsecolor"] = falsecolor_path
        attrs(file)["fname_truecolor"] = truecolor_path
        attrs(file)["iftversion"] = string(iftversion)
        attrs(file)["crs"] = latlondata[:crs]
        attrs(file)["crs_name"] = crs_name
        attrs(file)["reference"] = reference
        attrs(file)["contact"] = contact

        @info "Create group index"
        group_index = create_group(file, "index")
        group_index["time"] = ptsunix
        group_index["x"] = latlondata[:X]
        group_index["y"] = latlondata[:Y]

        @info "Create group floe_properties"
        group_floe_properties = create_group(file, "floe_properties")
        if nrow(props) > 0
            write_dataset(
                group_floe_properties, "properties", [copy(row) for row in eachrow(props)]
            )  # `copy(row)` converts the DataSetRow to a NamedTuple
            attrs(group_floe_properties)["Description of properties"] = """Area units (`area`, `convex_area`) are in sq. kilometers, length units (`minor_axis_length`, `major_axis_length`, and `perimeter`) in kilometers, and `orientation` in radians (see the description of properties attribute.) Latitude and longitude coordinates are in degrees, and the stereographic coordinates `x` and `y` are in meters relative to the $crs_name projection. """
        else
            attrs(group_floe_properties)["Description of properties"] = "No floes detected"
        end

        @info "Choose labeled data type"
        mx = maximum(labeled)
        T = choose_dtype(mx)

        @info "Write labeled image"
        labeled_rectified = T.(permutedims(labeled))
        label_data_obj, label_data_dtype = create_dataset(
            group_floe_properties, "labeled_image", labeled_rectified
        )
        attrs(label_data_obj)["CLASS"] = "IMAGE"
        attrs(label_data_obj)["IMAGE_SUBCLASS"] = "IMAGE_INDEXED"
        attrs(label_data_obj)["IMAGE_MINMAXRANGE"] = [
            minimum(labeled_rectified), maximum(labeled_rectified)
        ]
        attrs(label_data_obj)["description"] = "Connected components of the segmented floe image using a 3x3 structuring element. The property matrix consists of the properties of each connected component."
        write_dataset(label_data_obj, label_data_dtype, labeled_rectified)

        @info "Create group classifications"
        group_classifications = create_group(file, "classifications")

        @info "Write cloud mask"
        cloud_mask_rectified = T.(permutedims(cloud_mask))
        cloud_mask_obj, cloud_mask_dtype = create_dataset(
            group_classifications, "cloud_mask", cloud_mask_rectified
        )
        attrs(cloud_mask_obj)["CLASS"] = "IMAGE"
        attrs(cloud_mask_obj)["IMAGE_SUBCLASS"] = "IMAGE_GRAYSCALE"
        attrs(cloud_mask_obj)["IMAGE_MINMAXRANGE"] = [
            minimum(cloud_mask_rectified), maximum(cloud_mask_rectified)
        ]
        attrs(cloud_mask_obj)["description"] = "Cloud mask."
        write_dataset(cloud_mask_obj, cloud_mask_dtype, cloud_mask_rectified)

        @info "Write landmask"
        landmask_rectified = T.(permutedims(landmask))
        landmask_obj, landmask_dtype = create_dataset(
            group_classifications, "landmask", landmask_rectified
        )
        attrs(landmask_obj)["CLASS"] = "IMAGE"
        attrs(landmask_obj)["IMAGE_SUBCLASS"] = "IMAGE_GRAYSCALE"
        attrs(landmask_obj)["IMAGE_MINMAXRANGE"] = [
            minimum(landmask_rectified), maximum(landmask_rectified)
        ]
        attrs(landmask_obj)["description"] = "Land mask."
        write_dataset(landmask_obj, landmask_dtype, landmask_rectified)

        @info "Write coastal buffer mask"
        coastal_buffer_rectified = T.(permutedims(coastal_buffer_mask))
        coastal_buffer_obj, coastal_buffer_dtype = create_dataset(
            group_classifications, "coastal_buffer_mask", coastal_buffer_rectified
        )
        attrs(coastal_buffer_obj)["CLASS"] = "IMAGE"
        attrs(coastal_buffer_obj)["IMAGE_SUBCLASS"] = "IMAGE_GRAYSCALE"
        attrs(coastal_buffer_obj)["IMAGE_MINMAXRANGE"] = [
            minimum(coastal_buffer_rectified), maximum(coastal_buffer_rectified)
        ]
        attrs(coastal_buffer_obj)["description"] = "Coastal buffer mask. This mask is 1 for pixels within a specified distance of the coast, and 0 elsewhere."
        write_dataset(coastal_buffer_obj, coastal_buffer_dtype, coastal_buffer_rectified)

        @info "Write ice mask"
        ice_mask_rectified = T.(permutedims(ice_mask))
        ice_mask_obj, ice_mask_dtype = create_dataset(
            group_classifications, "ice_mask", ice_mask_rectified
        )
        attrs(ice_mask_obj)["CLASS"] = "IMAGE"
        attrs(ice_mask_obj)["IMAGE_SUBCLASS"] = "IMAGE_GRAYSCALE"
        attrs(ice_mask_obj)["IMAGE_MINMAXRANGE"] = [
            minimum(ice_mask_rectified), maximum(ice_mask_rectified)
        ]
        attrs(ice_mask_obj)["description"] = "Ice mask. This mask is 1 for pixels classified as ice, and 0 elsewhere."
        write_dataset(ice_mask_obj, ice_mask_dtype, ice_mask_rectified)
    end
end

end

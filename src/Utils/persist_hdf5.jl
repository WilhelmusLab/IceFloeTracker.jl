
module persist_hdf5
using IceFloeTracker
using HDF5, Images, Dates, TimeZones, DataFrames

export make_hdf5

function choose_dtype(mx::T) where {T<:Integer}
    types = [UInt8, UInt16, UInt32, UInt64, UInt128]
    for t_ in types
        if typemin(t_) <= mx <= typemax(t_)
            return t_
        end
    end
    throw("$mx cannot be represented by any of $types")
end

function make_hdf5(
    output_path::AbstractString,
    ;
    passtime::ZonedDateTime,
    crs_ref_image_path::AbstractString,
    truecolor_path::AbstractString,
    falsecolor_path::AbstractString,
    labels_map_path::AbstractString,
    cloud_mask_path::AbstractString,
    landmask_path::AbstractString,
    coastal_buffer_mask_path::AbstractString,
    iftversion::VersionNumber=pkgversion(IceFloeTracker),
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

    labeled = load(labels_map_path) |> channelview .|> x -> reinterpret(UInt16, x)
    props = regionprops_table(SegmentedImage(labeled, Int64.(labeled)))

    cloud_mask = load(cloud_mask_path) |> channelview .|> x -> reinterpret(UInt8, x)
    landmask = load(landmask_path) |> channelview .|> x -> reinterpret(UInt8, x)
    coastal_buffer_mask =
        load(coastal_buffer_mask_path) |> channelview .|> x -> reinterpret(UInt8, x)

    colstodrop = [:row_centroid, :col_centroid, :min_row, :min_col, :max_row, :max_col]
    converttounits!(props, latlondata, colstodrop)

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
        g = create_group(file, "index")
        g["time"] = ptsunix
        g["x"] = latlondata[:X]
        g["y"] = latlondata[:Y]

        @info "Create group floe_properties"
        g = create_group(file, "floe_properties")
        @show nrow(props)
        if nrow(props) > 0
            write_dataset(g, "properties", [copy(row) for row in eachrow(props)])  # `copy(row)` converts the DataSetRow to a NamedTuple
            attrs(g)["Description of properties"] = """Area units (`area`, `convex_area`) are in sq. kilometers, length units (`minor_axis_length`, `major_axis_length`, and `perimeter`) in kilometers, and `orientation` in radians (see the description of properties attribute.) Latitude and longitude coordinates are in degrees, and the stereographic coordinates `x` and `y` are in meters relative to the $crs_name projection. """
        else
            attrs(g)["Description of properties"] = "No floes detected"
        end

        @info "Choose labeled data type"
        mx = maximum(labeled)
        T = choose_dtype(mx)

        @info "Write labeled image"
        imgdata = T.(permutedims(labeled))
        obj, dtype = create_dataset(g, "labeled_image", imgdata)
        attrs(obj)["CLASS"] = "IMAGE"
        attrs(obj)["IMAGE_SUBCLASS"] = "IMAGE_INDEXED"
        attrs(obj)["IMAGE_MINMAXRANGE"] = [minimum(imgdata), maximum(imgdata)]

        @info "Add description"
        attrs(obj)["description"] = "Connected components of the segmented floe image using a 3x3 structuring element. The property matrix consists of the properties of each connected component."
        write_dataset(obj, dtype, imgdata)
    end
end

end
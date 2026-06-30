import HDF5: Group, File, create_dataset

"""
    IceFloeTracker.HDF5.V2(;
        passtime::ZonedDateTime,
        crs_ref_image_path::AbstractString,
        truecolor::AbstractMatrix{<:Union{RGB,RGBA}},
        falsecolor::AbstractMatrix{<:Union{RGB,RGBA}},
        labeled::AbstractMatrix,
        props::DataFrame,
        cloud_mask::AbstractMatrix,
        ice_mask::AbstractMatrix,
        landmask::AbstractMatrix,
        coastal_buffer_mask::AbstractMatrix,
        iftversion::VersionNumber = pkgversion(@__MODULE__),
        file_version::VersionNumber = VersionNumber("1.0.0"),
        reference::AbstractString = "https://doi.org/10.1016/j.rse.2019.111406",
        contact::AbstractString = "mmwilhelmus@brown.edu",
    )

An object with results from a single segmentation to be saved as an HDF5 file with [`save_hdf5`](@ref). 

Includes:

- References
  - `passtime`: the timepoint of the observation
  - `crs_ref_image_path`: the path to a georeferenced image
  - `truecolor`: the truecolor image
  - `falsecolor`: the falsecolor image
  - `iftversion`: the version of IceFloeTracker.jl used to save the file
  - `file_version`: the version of the file format (for this object, "1.0.0")
  - `reference`: a DOI for the dataset to which the file belongs
  - `contact`: contact information for the author
- Images
  - `labeled`: the labeled image of connected components
  - `cloud_mask`: the cloud mask
  - `ice_mask`: the ice mask
  - `landmask`: the land mask
  - `coastal_buffer_mask`: the coastal buffer mask
- DataFrames
  - `props`: the measured properties of the floes

"""
@kwdef struct V2
    passtime::ZonedDateTime
    crs_ref_image_path::AbstractString
    truecolor::AbstractMatrix{<:Union{RGB,RGBA}}
    falsecolor::AbstractMatrix{<:Union{RGB,RGBA}}
    labeled::AbstractMatrix
    props::DataFrame
    cloud_mask::AbstractMatrix
    ice_mask::AbstractMatrix
    landmask::AbstractMatrix
    coastal_buffer_mask::AbstractMatrix
    iftversion::VersionNumber = pkgversion(@__MODULE__)
    file_version::VersionNumber = VersionNumber("2.0.0")
    reference::AbstractString = "https://doi.org/10.1016/j.rse.2019.111406"
    contact::AbstractString = "mmwilhelmus@brown.edu"
end

"""
    save_hdf5(path, V2(args...))

Write the [`V2`](@ref) object to storage.

The structure is:
🗂️ HDF5.File:
├─ 🏷️ file_version
├─ 🏷️ iftversion
├─ 🏷️ contact
├─ 🏷️ reference
├─ 🏷️ crs
├─ 🏷️ crs_name
├─ 🏷️ fname_falsecolor
├─ 🏷️ fname_truecolor
├─ 📂 classifications
│  ├─ 🔢 cloud_mask
│  │  ├─ 🏷️ CLASS
│  │  ├─ 🏷️ IMAGE_MINMAXRANGE
│  │  ├─ 🏷️ IMAGE_SUBCLASS
│  │  └─ 🏷️ description
│  ├─ 🔢 coastal_buffer_mask
│  │  ├─ 🏷️ CLASS
│  │  ├─ 🏷️ IMAGE_MINMAXRANGE
│  │  ├─ 🏷️ IMAGE_SUBCLASS
│  │  └─ 🏷️ description
│  ├─ 🔢 ice_mask
│  │  ├─ 🏷️ CLASS
│  │  ├─ 🏷️ IMAGE_MINMAXRANGE
│  │  ├─ 🏷️ IMAGE_SUBCLASS
│  │  └─ 🏷️ description
│  └─ 🔢 landmask
│     ├─ 🏷️ CLASS
│     ├─ 🏷️ IMAGE_MINMAXRANGE
│     ├─ 🏷️ IMAGE_SUBCLASS
│     └─ 🏷️ description
├─ 📂 floe_properties
│  ├─ 🏷️ Description of properties
│  ├─ 🔢 labeled_image
│  │  ├─ 🏷️ CLASS
│  │  ├─ 🏷️ IMAGE_MINMAXRANGE
│  │  ├─ 🏷️ IMAGE_SUBCLASS
│  │  └─ 🏷️ description
│  └─ 🔢 properties
└─ 📂 index
   ├─ 🔢 time
   ├─ 🔢 x
   └─ 🔢 y
"""
function save_hdf5(output_path::AbstractString, s::V2;)
    ptsunix = Int64(Dates.datetime2unix(DateTime(s.passtime)))
    latlondata = latlon(s.crs_ref_image_path)

    crs_code = latlondata[:crs]
    crs_name = get_crs_name(crs_code)
    projection_dataset_name = get_projection_name(crs_code)

    h5open(output_path, "w") do file
        attrs(file)["file_version"] = string(s.file_version)
        attrs(file)["iftversion"] = string(s.iftversion)
        attrs(file)["reference"] = s.reference
        attrs(file)["contact"] = s.contact

        dset = create_dataset(file, projection_dataset_name, String, (1,))
        attrs(dset)["name"] = crs_name
        attrs(dset)["crs_wkt"] = latlondata[:crs_wkt]
        attrs(dset)["spatial_ref"] = latlondata[:crs_wkt]
        attrs(dset)["long_name"] = "CRS Definition"
        attrs(dset)["GeoTransform"] = join(Int64.(latlondata[:geotransform]), " ")

        file["x"] = latlondata[:X]
        attrs(file["x"])["standard_name"] = "projection_x_coordinate"
        attrs(file["x"])["long_name"] = "x coordinate of projection"
        attrs(file["x"])["units"] = "m"

        file["y"] = latlondata[:Y]
        attrs(file["y"])["standard_name"] = "projection_y_coordinate"
        attrs(file["y"])["long_name"] = "y coordinate of projection"
        attrs(file["y"])["units"] = "m"

        file["time"] = ptsunix

        create_color_dataset(
            file, "falsecolor", s.falsecolor, "Falsecolor image", projection_dataset_name
        )
        create_color_dataset(
            file, "truecolor", s.truecolor, "Truecolor image", projection_dataset_name
        )

        floe_properties_group = create_group(file, "floe_properties")
        create_floe_properties_dataset(
            floe_properties_group, s.props, "floe-properties", crs_name
        )

        create_labeled_dataset(
            file,
            "labeled_image",
            s.labeled,
            "Connected components of the segmented floe image using a 3x3 structuring element. The property matrix consists of the properties of each connected component.",
            projection_dataset_name,
        )

        group_classifications = create_group(file, "classifications")

        create_mask_dataset(
            group_classifications,
            "cloud_mask",
            s.cloud_mask,
            "Cloud mask. This mask is 1 for pixels classified as cloud, and 0 elsewhere.",
            projection_dataset_name,
        )

        create_mask_dataset(
            group_classifications,
            "landmask",
            s.landmask,
            "Land mask. This mask is 1 for pixels classified as land, and 0 elsewhere.",
            projection_dataset_name,
        )

        create_mask_dataset(
            group_classifications,
            "coastal_buffer_mask",
            s.coastal_buffer_mask,
            "Coastal buffer mask. This mask is 1 for pixels within a specified distance of the coast, and 0 elsewhere.",
            projection_dataset_name,
        )

        create_mask_dataset(
            group_classifications,
            "ice_mask",
            s.ice_mask,
            "Ice mask. This mask is 1 for pixels classified as ice, and 0 elsewhere.",
            projection_dataset_name,
        )
    end
    return nothing
end

function _load_v2(file)
    passtime = ZonedDateTime(unix2datetime(read(file["index/time"])), tz"UTC")
    crs_ref_image_path = attrs(file)["fname_truecolor"]
    truecolor_path = attrs(file)["fname_truecolor"]
    falsecolor_path = attrs(file)["fname_falsecolor"]
    labeled = read(file["floe_properties/labeled_image"]) .|> Int
    props = DataFrame(read(file["floe_properties/properties"]))
    landmask = read(file["classifications/landmask"]) |> binarize_mask .|> Gray
    cloud_mask = read(file["classifications/cloud_mask"]) |> binarize_mask .|> Gray
    ice_mask = read(file["classifications/ice_mask"]) |> binarize_mask .|> Gray
    coastal_buffer_mask =
        read(file["classifications/coastal_buffer_mask"]) |> binarize_mask .|> Gray
    iftversion = VersionNumber(attrs(file)["iftversion"])
    reference = attrs(file)["reference"]
    contact = attrs(file)["contact"]

    return V2(;
        passtime,
        crs_ref_image_path,
        truecolor_path,
        falsecolor_path,
        labeled,
        props,
        cloud_mask,
        ice_mask,
        landmask,
        coastal_buffer_mask,
        iftversion,
        reference,
        contact,
    )
end

function get_crs_name(crs_code)
    crs_dict = Dict(
        3413 => "EPSG:3413 NSIDC north polar stereographic",
        3031 => "EPSG:3031 NSIDC south polar stereographic",
        4326 => "EPSG:4326 WGS84 lat/lon",
        3857 => "EPSG:3857 Web Mercator",
    )
    crs_name = get(crs_dict, crs_code) do
        crs_name_ = "EPSG:$(string(crs_code))"
        @warn "CRS $crs_code not recognized. CRS will be recorded as $crs_name_ in the output file attributes, but no short name will be provided."
        return crs_name_
    end
    return crs_name
end

function get_projection_name(crs_code)::String
    @show crs_code
    crs_short_name_dict = Dict(
        3413 => "north_polar_stereographic",
        3031 => "south_polar_stereographic",
        4326 => "wgs84_lat_lon",
        3857 => "web_mercator",
    )
    projection_dataset_name = get(crs_short_name_dict, crs_code) do
        crs_name_ = "geolocation"
        return crs_name_
    end
    return projection_dataset_name
end

function create_floe_properties_dataset(
    group::Union{File,Group},
    props::DataFrame,
    name::AbstractString="floe-properties",
    crs_name::AbstractString="",
)
    @show props
    props_ = convert_missing_to_nan(props)
    @show props_
    if nrow(props_) > 0
        create_dataset(group, name, [copy(row) for row in eachrow(props_)])  # `copy(row)` converts the DataSetRow to a NamedTuple
        attrs(group)["Description"] = """Area units (`area`, `convex_area`) are in sq. kilometers, length units (`minor_axis_length`, `major_axis_length`, and `perimeter`) in kilometers, and `orientation` in radians (see the description of properties attribute.) Latitude and longitude coordinates are in degrees, and the stereographic coordinates `x` and `y` are in meters relative to the $crs_name projection. """
    else
        create_dataset(group, name, [])
        attrs(group)["Description"] = "No floes detected"
    end
end

function create_mask_dataset(
    group::Union{File,Group},
    name::AbstractString,
    mask::AbstractMatrix,
    description::AbstractString="",
    projection_dataset_name::AbstractString="geolocation",
)
    mx = maximum(mask)
    T = choose_dtype(mx)
    mask_rectified = T.(permutedims(mask, (2, 1)))
    mask_obj, mask_dtype = create_dataset(group, name, mask_rectified)
    attrs(mask_obj)["CLASS"] = "IMAGE"
    attrs(mask_obj)["IMAGE_SUBCLASS"] = "IMAGE_GRAYSCALE"
    attrs(mask_obj)["IMAGE_MINMAXRANGE"] = [
        minimum(mask_rectified), maximum(mask_rectified)
    ]
    attrs(mask_obj)["description"] = description
    attrs(mask_obj)["grid_mapping"] = projection_dataset_name
    write_dataset(mask_obj, mask_dtype, mask_rectified)
end

function create_labeled_dataset(
    group::Union{File,Group},
    name::AbstractString,
    labeled::AbstractMatrix,
    description::AbstractString="",
    projection_dataset_name::AbstractString="geolocation",
)
    mx = maximum(labeled)
    T = choose_dtype(mx)
    labeled_rectified = T.(permutedims(labeled, (2, 1)))
    label_data_obj, label_data_dtype = create_dataset(group, name, labeled_rectified)
    attrs(label_data_obj)["CLASS"] = "IMAGE"
    attrs(label_data_obj)["IMAGE_SUBCLASS"] = "IMAGE_INDEXED"
    attrs(label_data_obj)["IMAGE_MINMAXRANGE"] = [
        minimum(labeled_rectified), maximum(labeled_rectified)
    ]
    attrs(label_data_obj)["description"] = description
    attrs(label_data_obj)["grid_mapping"] = projection_dataset_name
    write_dataset(label_data_obj, label_data_dtype, labeled_rectified)
end

function create_color_dataset(
    group::Union{File,Group},
    name::AbstractString,
    img::AbstractMatrix{<:Union{RGB,RGBA}},
    description::AbstractString="",
    projection_dataset_name::AbstractString="geolocation",
)
    img_rectified = permutedims(rawview(channelview(img)), (3, 2, 1))
    el = eltype(img_rectified)
    img_obj, img_dtype = create_dataset(group, name, img_rectified)
    attrs(img_obj)["CLASS"] = "IMAGE"
    attrs(img_obj)["IMAGE_SUBCLASS"] = "IMAGE_TRUECOLOR"
    attrs(img_obj)["IMAGE_VERSION"] = "1.2"
    attrs(img_obj)["INTERLACE_MODE"] = "INTERLACE_PLANE"
    attrs(img_obj)["IMAGE_MINMAXRANGE"] = [el(typemin(el)), el(typemax(el))]
    attrs(img_obj)["description"] = description
    attrs(img_obj)["grid_mapping"] = projection_dataset_name

    write_dataset(img_obj, img_dtype, img_rectified)
end
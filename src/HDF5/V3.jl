"""
    IceFloeTracker.HDF5.V3(;
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
        file_version::VersionNumber = VersionNumber("3.0.0"),
        reference::AbstractString = "https://doi.org/10.1016/j.rse.2019.111406",
        contact::AbstractString = "mmwilhelmus@brown.edu",
    )

An object with results from a single segmentation to be saved as a proper
netCDF-4 file with [`save_hdf5`](@ref).

Includes:

- References
  - `passtime`: the timepoint of the observation
  - `crs_ref_image_path`: the path to a georeferenced image
  - `truecolor`: the truecolor image
  - `falsecolor`: the falsecolor image
  - `iftversion`: the version of IceFloeTracker.jl used to save the file
  - `file_version`: the version of the file format (for this object, "3.0.0")
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
@kwdef struct V3
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
    file_version::VersionNumber = VersionNumber("3.0.0")
    reference::AbstractString = "https://doi.org/10.1016/j.rse.2019.111406"
    contact::AbstractString = "mmwilhelmus@brown.edu"
end

"""
    save_hdf5(path, V3(args...))

Write the [`V3`](@ref) object to storage as a proper netCDF-4 file.

The file uses `NCDatasets` and follows CF conventions. Image variables carry
`CLASS`/`IMAGE_SUBCLASS` attributes so viewers that recognise those conventions
can display them directly. The `x`, `y`, and `time` dimensions are defined at
the root level and inherited by all groups.

Structure:

```
📦 netCDF-4 file
├─ 🏷️ file_version, iftversion, contact, reference   (global attributes)
├─ 🔢 geolocation   (scalar Int32, CRS grid-mapping variable)
│  └─ 🏷️ name, crs_wkt, spatial_ref, long_name, GeoTransform
├─ 🔢 x(x)          (projection x / easting coordinates, metres)
├─ 🔢 y(y)          (projection y / northing coordinates, metres)
├─ 🔢 time(time)    (seconds since 1970-01-01)
├─ 📂 inputs
│  ├─ 🔢 falsecolor(band, y, x)
│  │  └─ 🏷️ CLASS, IMAGE_MINMAXRANGE, IMAGE_SUBCLASS, IMAGE_VERSION,
│  │        INTERLACE_MODE, description, grid_mapping
│  └─ 🔢 truecolor(band, y, x)
│     └─ (same attributes as falsecolor)
├─ 📂 classifications
│  ├─ 🔢 labeled_image(y, x)
│  │  └─ 🏷️ CLASS, IMAGE_MINMAXRANGE, IMAGE_SUBCLASS, description, grid_mapping
│  ├─ 🔢 cloud_mask(y, x)
│  ├─ 🔢 coastal_buffer_mask(y, x)
│  ├─ 🔢 ice_mask(y, x)
│  └─ 🔢 landmask(y, x)
└─ 📂 floe-properties
   ├─ 🏷️ Description, label_variable, labeled_image   (group attributes)
   ├─ 🔢 label(floe)   (links pixel values in labeled_image to rows here)
   └─ 🔢 <column>(floe) for every other column in the props DataFrame
```

"""
function save_hdf5(output_path::AbstractString, s::V3;)
    ptsunix = Int64(Dates.datetime2unix(DateTime(s.passtime)))
    latlondata = latlon(s.crs_ref_image_path)

    crs_code = latlondata[:crs]
    crs_name = get_crs_name(crs_code)
    projection_dataset_name = "geolocation"

    X = Float64.(latlondata[:X])
    Y = Float64.(latlondata[:Y])
    nx = length(X)
    ny = length(Y)

    NCDataset(output_path, "c") do ds
        # Global attributes
        ds.attrib["file_version"] = string(s.file_version)
        ds.attrib["iftversion"] = string(s.iftversion)
        ds.attrib["reference"] = s.reference
        ds.attrib["contact"] = s.contact

        # Dimensions defined at root are inherited by all groups
        defDim(ds, "x", nx)
        defDim(ds, "y", ny)
        defDim(ds, "time", 1)

        # CRS / geolocation grid_mapping variable (scalar)
        vcrs = defVar(ds, projection_dataset_name, Int32, ())
        vcrs.attrib["name"] = crs_name
        vcrs.attrib["crs_wkt"] = latlondata[:crs_wkt]
        vcrs.attrib["spatial_ref"] = latlondata[:crs_wkt]
        vcrs.attrib["long_name"] = "CRS Definition"
        vcrs.attrib["GeoTransform"] = join(Int64.(latlondata[:geotransform]), " ")
        vcrs[] = Int32(0)

        # x coordinate variable
        vx = defVar(ds, "x", Float64, ("x",))
        vx.attrib["standard_name"] = "projection_x_coordinate"
        vx.attrib["long_name"] = "x coordinate of projection"
        vx.attrib["units"] = "m"
        vx[:] = X

        # y coordinate variable
        vy = defVar(ds, "y", Float64, ("y",))
        vy.attrib["standard_name"] = "projection_y_coordinate"
        vy.attrib["long_name"] = "y coordinate of projection"
        vy.attrib["units"] = "m"
        vy[:] = Y

        # time variable
        vt = defVar(ds, "time", Int64, ("time",))
        vt.attrib["units"] = "seconds since 1970-01-01 00:00:00"
        vt.attrib["calendar"] = "gregorian"
        vt.attrib["long_name"] = "time of satellite overpass"
        vt[:] = [ptsunix]

        # inputs group: colour imagery
        grp_inputs = defGroup(ds, "inputs")
        nchannels = size(channelview(s.truecolor), 1)
        defDim(grp_inputs, "band", nchannels)
        nc_create_color_dataset(
            grp_inputs,
            "falsecolor",
            s.falsecolor,
            "Falsecolor image",
            projection_dataset_name,
        )
        nc_create_color_dataset(
            grp_inputs, "truecolor", s.truecolor, "Truecolor image", projection_dataset_name
        )

        # classifications group: masks + labelled segmentation map
        grp_class = defGroup(ds, "classifications")
        nc_create_labeled_dataset(
            grp_class,
            "labeled_image",
            s.labeled,
            "Connected components of the segmented floe image using a 3x3 " *
            "structuring element. The property matrix consists of the properties " *
            "of each connected component.",
            projection_dataset_name,
        )
        nc_create_mask_dataset(
            grp_class,
            "cloud_mask",
            s.cloud_mask,
            "Cloud mask. This mask is 1 for pixels classified as cloud, and 0 elsewhere.",
            projection_dataset_name,
        )
        nc_create_mask_dataset(
            grp_class,
            "landmask",
            s.landmask,
            "Land mask. This mask is 1 for pixels classified as land, and 0 elsewhere.",
            projection_dataset_name,
        )
        nc_create_mask_dataset(
            grp_class,
            "coastal_buffer_mask",
            s.coastal_buffer_mask,
            "Coastal buffer mask. This mask is 1 for pixels within a specified " *
            "distance of the coast, and 0 elsewhere.",
            projection_dataset_name,
        )
        nc_create_mask_dataset(
            grp_class,
            "ice_mask",
            s.ice_mask,
            "Ice mask. This mask is 1 for pixels classified as ice, and 0 elsewhere.",
            projection_dataset_name,
        )

        # floe-properties group: one variable per property column
        grp_props = defGroup(ds, "floe-properties")
        nc_create_floe_properties(
            grp_props, s.props, crs_name, "classifications/labeled_image"
        )
    end

    # The netCDF-C library reserves the "CLASS" attribute name (it is used
    # internally by HDF5 for typed datasets).  We set it after closing the
    # NCDatasets handle by writing directly to the underlying HDF5 file.
    image_paths = [
        "/inputs/falsecolor",
        "/inputs/truecolor",
        "/classifications/labeled_image",
        "/classifications/cloud_mask",
        "/classifications/landmask",
        "/classifications/coastal_buffer_mask",
        "/classifications/ice_mask",
    ]
    h5open(output_path, "r+") do file
        for path in image_paths
            attrs(file[path])["CLASS"] = "IMAGE"
        end
        API.h5l_create_hard(
            file.id,
            "/classifications/labeled_image",
            file.id,
            "/floe-properties/labeled_image",
            API.H5P_DEFAULT,
            API.H5P_DEFAULT,
        )
    end

    return nothing
end

"""
    nc_create_mask_dataset(grp, name, mask, description, projection_dataset_name)

Define a `UInt8` binary mask variable with dimensions `(y, x)` in the NCDatasets
group `grp`. The parent group must supply the `x` and `y` dimensions. Standard
image attributes (`CLASS`, `IMAGE_SUBCLASS`, `IMAGE_MINMAXRANGE`, `grid_mapping`,
`description`) are attached to the variable.
"""
function nc_create_mask_dataset(
    grp::NCDataset,
    name::AbstractString,
    mask::AbstractMatrix,
    description::AbstractString="",
    projection_dataset_name::AbstractString="geolocation",
)
    mask_uint8 = permutedims(UInt8.(Bool.(mask)), (2, 1))  # (nx, ny)
    v = defVar(grp, name, UInt8, ("x", "y"))
    v.attrib["IMAGE_SUBCLASS"] = "IMAGE_GRAYSCALE"
    v.attrib["IMAGE_MINMAXRANGE"] = UInt8[minimum(mask_uint8), maximum(mask_uint8)]
    v.attrib["description"] = description
    v.attrib["grid_mapping"] = projection_dataset_name
    v[:, :] = mask_uint8
    return v
end

"""
    nc_create_labeled_dataset(grp, name, labeled, description, projection_dataset_name)

Define an integer-typed segmentation label variable with dimensions `(y, x)` in
the NCDatasets group `grp`. The integer type is the smallest that can represent
the maximum label value. The parent group must supply the `x` and `y` dimensions.
"""
function nc_create_labeled_dataset(
    grp::NCDataset,
    name::AbstractString,
    labeled::AbstractMatrix,
    description::AbstractString="",
    projection_dataset_name::AbstractString="geolocation",
)
    mx = Int64(maximum(labeled))
    T = choose_dtype(mx)
    labeled_T = permutedims(T.(labeled), (2, 1))  # (nx, ny)
    v = defVar(grp, name, T, ("x", "y"))
    v.attrib["IMAGE_SUBCLASS"] = "IMAGE_INDEXED"
    v.attrib["IMAGE_MINMAXRANGE"] = T[T(minimum(labeled)), maximum(labeled_T)]
    v.attrib["description"] = description
    v.attrib["grid_mapping"] = projection_dataset_name
    v[:, :] = labeled_T
    return v
end

"""
    nc_create_color_dataset(grp, name, img, description, projection_dataset_name)

Define a `UInt8` colour image variable with dimensions `(band, y, x)` in the
NCDatasets group `grp`. The group must already define a `band` dimension, and the
parent group must supply the `x` and `y` dimensions.
"""
function nc_create_color_dataset(
    grp::NCDataset,
    name::AbstractString,
    img::AbstractMatrix{<:Union{RGB,RGBA}},
    description::AbstractString="",
    projection_dataset_name::AbstractString="geolocation",
)
    img_raw = permutedims(UInt8.(rawview(channelview(img))), (3, 2, 1))  # (nx, ny, nchannels)
    v = defVar(grp, name, UInt8, ("x", "y", "band"))
    v.attrib["IMAGE_SUBCLASS"] = "IMAGE_TRUECOLOR"
    v.attrib["IMAGE_VERSION"] = "1.2"
    v.attrib["INTERLACE_MODE"] = "INTERLACE_PLANE"
    v.attrib["IMAGE_MINMAXRANGE"] = UInt8[0, 255]
    v.attrib["description"] = description
    v.attrib["grid_mapping"] = projection_dataset_name
    v[:, :, :] = img_raw
    return v
end

"""
    nc_create_floe_properties(grp, props, crs_name, labeled_image_path)

Define one netCDF variable per column in the `props` DataFrame, all sharing a
`floe` dimension. Integer columns are stored as `Int64`; floating-point columns
as `Float64` with a `NaN` fill value; other columns as `String`. The `label`
variable links each floe back to pixel values in the labeled image via the group
attributes `label_variable` and `labeled_image`.
"""
function nc_create_floe_properties(
    grp::NCDataset,
    props::DataFrame,
    crs_name::AbstractString="",
    labeled_image_path::AbstractString="classifications/labeled_image",
)
    props_ = convert_missing_to_nan(props)
    nfloes = nrow(props_)

    defDim(grp, "floe", nfloes)

    grp.attrib["Description"] = (
        "Area units (`area`, `convex_area`) are in sq. kilometers, " *
        "length units (`minor_axis_length`, `major_axis_length`, and `perimeter`) " *
        "in kilometers, and `orientation` in radians. " *
        "Latitude and longitude coordinates are in degrees, and the " *
        "stereographic coordinates `x` and `y` are in metres relative to the " *
        "$crs_name projection."
    )
    grp.attrib["label_variable"] = "label"
    grp.attrib["labeled_image"] = labeled_image_path

    for col_name in names(props_)
        col_data = props_[!, col_name]
        T = eltype(col_data)

        v = if T <: Integer
            defVar(grp, col_name, Int64, ("floe",))
        elseif T <: AbstractFloat
            defVar(grp, col_name, Float64, ("floe",); fillvalue=NaN)
        else
            defVar(grp, col_name, String, ("floe",))
        end

        if col_name == "label"
            v.attrib["long_name"] = "floe label"
            v.attrib["comment"] = "Pixel values in $labeled_image_path equal these labels"
        end

        if nfloes > 0
            v[:] = col_data
        end
    end
    return nothing
end
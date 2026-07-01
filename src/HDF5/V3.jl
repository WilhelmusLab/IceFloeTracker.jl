"""
    IceFloeTracker.HDF5.V3(;
        passtime::ZonedDateTime,
        crs::NamedTuple,
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
  - `crs`: the CRS and geospatial data for the image, as returned by [`latlon`](@ref)
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
    crs::NamedTuple
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
├─ 🏷️ Description, label_variable, labeled_image     (floe-properties attributes)
├─ 🔢 geolocation   (scalar Int32, CRS grid-mapping variable)
│  └─ 🏷️ name, crs_wkt, spatial_ref, long_name, GeoTransform
├─ 🔢 x(x)          (projection x / easting coordinates, metres)
├─ 🔢 y(y)          (projection y / northing coordinates, metres)
├─ 🔢 time(time)    (seconds since 1970-01-01)
├─ 🔢 falsecolor(x, y, band_falsecolor)
│  └─ 🏷️ CLASS, IMAGE_MINMAXRANGE, IMAGE_SUBCLASS, IMAGE_VERSION,
│        INTERLACE_MODE, description, grid_mapping
├─ 🔢 truecolor(x, y, band_truecolor)
│  └─ (same attributes as falsecolor)
├─ 🔢 labeled_image(x, y)
│  └─ 🏷️ CLASS, IMAGE_MINMAXRANGE, IMAGE_SUBCLASS, description, grid_mapping
├─ 🔢 cloud_mask(x, y)
├─ 🔢 coastal_buffer_mask(x, y)
├─ 🔢 ice_mask(x, y)
├─ 🔢 landmask(x, y)
├─ 🔢 label(floe)   (links pixel values in labeled_image to rows here)
└─ 🔢 <column>(floe) for every other column in the props DataFrame
```

"""
function save_hdf5(output_path::AbstractString, s::V3;)
    ptsunix = Int64(Dates.datetime2unix(DateTime(s.passtime)))
    latlondata = s.crs

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
        vcrs.attrib["EPSG"] = Int32(latlondata[:crs])
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

        # colour imagery
        nchannels_tc = size(channelview(s.truecolor), 1)
        nchannels_fc = size(channelview(s.falsecolor), 1)
        defDim(ds, "band_truecolor", nchannels_tc)
        defDim(ds, "band_falsecolor", nchannels_fc)
        nc_create_color_dataset(
            ds,
            "falsecolor",
            s.falsecolor,
            "Falsecolor image",
            projection_dataset_name,
            "band_falsecolor",
        )
        nc_create_color_dataset(
            ds,
            "truecolor",
            s.truecolor,
            "Truecolor image",
            projection_dataset_name,
            "band_truecolor",
        )

        # masks + labelled segmentation map
        nc_create_labeled_dataset(
            ds,
            "labeled_image",
            s.labeled,
            "Connected components of the segmented floe image using a 3x3 " *
            "structuring element. The property matrix consists of the properties " *
            "of each connected component.",
            projection_dataset_name,
        )
        nc_create_mask_dataset(
            ds,
            "cloud_mask",
            s.cloud_mask,
            "Cloud mask. This mask is 1 for pixels classified as cloud, and 0 elsewhere.",
            projection_dataset_name,
        )
        nc_create_mask_dataset(
            ds,
            "landmask",
            s.landmask,
            "Land mask. This mask is 1 for pixels classified as land, and 0 elsewhere.",
            projection_dataset_name,
        )
        nc_create_mask_dataset(
            ds,
            "coastal_buffer_mask",
            s.coastal_buffer_mask,
            "Coastal buffer mask. This mask is 1 for pixels within a specified " *
            "distance of the coast, and 0 elsewhere.",
            projection_dataset_name,
        )
        nc_create_mask_dataset(
            ds,
            "ice_mask",
            s.ice_mask,
            "Ice mask. This mask is 1 for pixels classified as ice, and 0 elsewhere.",
            projection_dataset_name,
        )

        # floe properties: one variable per property column
        nc_create_floe_properties(ds, s.props, "labeled_image")
    end

    # The netCDF-C library reserves the "CLASS" attribute name (it is used
    # internally by HDF5 for typed datasets).  We set it after closing the
    # NCDatasets handle by writing directly to the underlying HDF5 file.
    image_paths = [
        "/falsecolor",
        "/truecolor",
        "/labeled_image",
        "/cloud_mask",
        "/landmask",
        "/coastal_buffer_mask",
        "/ice_mask",
    ]
    h5open(output_path, "r+") do file
        for path in image_paths
            attrs(file[path])["CLASS"] = "IMAGE"
        end
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
    nc_create_color_dataset(grp, name, img, description, projection_dataset_name, band_dim_name)

Define a `UInt8` colour image variable with dimensions `(x, y, band_dim_name)` in the
NCDatasets group `grp`. The group must already define a dimension named `band_dim_name`,
and the parent group must supply the `x` and `y` dimensions.
"""
function nc_create_color_dataset(
    grp::NCDataset,
    name::AbstractString,
    img::AbstractMatrix{<:Union{RGB,RGBA}},
    description::AbstractString="",
    projection_dataset_name::AbstractString="geolocation",
    band_dim_name::AbstractString="band",
)
    img_raw = permutedims(UInt8.(rawview(channelview(img))), (3, 2, 1))  # (nx, ny, nchannels)
    v = defVar(grp, name, UInt8, ("x", "y", band_dim_name))
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
attributes `label_variable` and `labeled_image`. CF `units`, `long_name`, and
`standard_name` attributes are attached to each recognized column.
"""
function nc_create_floe_properties(
    grp::NCDataset, props::DataFrame, labeled_image_path::AbstractString="labeled_image"
)
    props_ = convert_missing_to_nan(props)
    nfloes = nrow(props_)

    defDim(grp, "floe", nfloes)

    grp.attrib["label_variable"] = "label"
    grp.attrib["labeled_image"] = labeled_image_path

    # "x" and "y" are already taken by the root-level coordinate variables;
    # remap the floe stereographic-coordinate columns to avoid the name clash.
    col_name_map = Dict("x" => "x_floe", "y" => "y_floe")

    # CF units and metadata keyed by the original DataFrame column name.
    col_attrs = Dict(
        "label" => (
            long_name="floe label",
            comment="Pixel values in $labeled_image_path equal these labels",
        ),
        "area" => (units="km^2", long_name="floe area"),
        "convex_area" => (units="km^2", long_name="convex hull area"),
        "minor_axis_length" => (units="km", long_name="minor axis length"),
        "major_axis_length" => (units="km", long_name="major axis length"),
        "perimeter" => (units="km", long_name="perimeter"),
        "orientation" => (units="rad", long_name="orientation angle"),
        "latitude" => (
            units="degrees_north",
            standard_name="latitude",
            long_name="latitude of floe centroid",
        ),
        "longitude" => (
            units="degrees_east",
            standard_name="longitude",
            long_name="longitude of floe centroid",
        ),
        "x" => (
            units="m",
            standard_name="projection_x_coordinate",
            long_name="x coordinate of floe centroid",
        ),
        "y" => (
            units="m",
            standard_name="projection_y_coordinate",
            long_name="y coordinate of floe centroid",
        ),
    )

    for col_name in names(props_)
        nc_name = get(col_name_map, col_name, col_name)
        col_data = props_[!, col_name]
        T = eltype(col_data)

        v = if T <: Integer
            defVar(grp, nc_name, Int64, ("floe",))
        elseif T <: AbstractFloat
            defVar(grp, nc_name, Float64, ("floe",); fillvalue=NaN)
        else
            defVar(grp, nc_name, String, ("floe",))
        end

        if haskey(col_attrs, col_name)
            for (k, val) in pairs(col_attrs[col_name])
                v.attrib[string(k)] = val
            end
        end

        if nfloes > 0
            v[:] = col_data
        end
    end
    return nothing
end

"""
    _load_v3(input_path)

Load a V3 netCDF-4 file written by [`save_hdf5`](@ref) and return a [`V3`](@ref) object.
"""
function _load_v3(input_path::AbstractString)
    NCDataset(input_path, "r") do ds
        # Global attributes
        iftversion = VersionNumber(ds.attrib["iftversion"])
        file_version = VersionNumber(ds.attrib["file_version"])
        reference = ds.attrib["reference"]
        contact = ds.attrib["contact"]

        # Reconstruct CRS data from the geolocation variable
        geoloc = ds["geolocation"]
        crs_wkt = geoloc.attrib["crs_wkt"]
        crs_code = Int64(geoloc.attrib["EPSG"])
        geotransform = parse.(Float64, split(geoloc.attrib["GeoTransform"]))
        X = Array(ds["x"].var)
        Y = Array(ds["y"].var)
        crs = (crs=crs_code, crs_wkt=crs_wkt, X=X, Y=Y, geotransform=geotransform)

        # passtime — raw Int64 stored as "seconds since 1970-01-01"
        ptsunix = Int64(ds["time"].var[1])
        passtime = ZonedDateTime(Dates.unix2datetime(ptsunix), tz"UTC")

        # Colour images: stored as (nx, ny, nb) UInt8 → reconstruct Matrix{RGB/RGBA}
        function read_color(varname)
            raw = collect(permutedims(Array(ds[varname].var), (3, 2, 1)))  # (nb, ny, nx)
            nb = size(raw, 1)
            n0f8 = reinterpret(N0f8, raw)
            return if nb == 3
                collect(colorview(RGB{N0f8}, n0f8))
            else
                collect(colorview(RGBA{N0f8}, n0f8))
            end
        end
        truecolor = read_color("truecolor")
        falsecolor = read_color("falsecolor")

        # Labeled image: stored as (nx, ny) → (ny, nx)
        labeled = permutedims(Array(ds["labeled_image"].var), (2, 1)) .|> Int

        # Masks: stored as (nx, ny) UInt8 → (ny, nx) BitMatrix
        function read_mask(varname)
            return permutedims(Array(ds[varname].var), (2, 1)) |> binarize_mask
        end
        cloud_mask = read_mask("cloud_mask")
        landmask = read_mask("landmask")
        ice_mask = read_mask("ice_mask")
        coastal_buffer_mask = read_mask("coastal_buffer_mask")

        # Floe properties: all variables sharing the "floe" dimension, in file order.
        # Reverse the x_floe/y_floe name remapping applied during save.
        col_name_remap = Dict("x_floe" => "x", "y_floe" => "y")
        prop_cols = Pair{String,AbstractVector}[]
        for varname in keys(ds)
            v = ds[varname]
            "floe" ∈ dimnames(v) || continue
            col_name = get(col_name_remap, varname, varname)
            # Use .var to bypass fill-value masking so NaN is preserved as NaN
            push!(prop_cols, col_name => Array(v.var))
        end
        props = isempty(prop_cols) ? DataFrame() : DataFrame(prop_cols)

        return V3(;
            passtime,
            crs,
            truecolor,
            falsecolor,
            labeled,
            props,
            cloud_mask,
            ice_mask,
            landmask,
            coastal_buffer_mask,
            iftversion,
            file_version,
            reference,
            contact,
        )
    end
end
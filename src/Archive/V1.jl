using ArchGDAL: importWKT, importEPSG

TRUECOLOR_BANDLABELS = [
    "R=band_1 (red, 0.620–0.670 µm)",
    "G=band_4 (green, 0.545–0.565 µm)",
    "B=band_3 (blue, 0.459–0.479 µm)",
]

FALSECOLOR_BANDLABELS = [
    "R=band_7 (mid-IR, 2.105–2.155 µm)",
    "G=band_2 (NIR, 0.841–0.876 µm)",
    "B=band_1 (red, 0.620–0.670 µm)",
]

COLUMN_ATTRIBUTES = Dict(
    "label" => (
        long_name="floe label",
        comment="Pixel values in labeled_image equal these labels",
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

"""
    IceFloeTracker.Archive.V1(;
        passtime::ZonedDateTime,
        crs::NamedTuple,
        modis_truecolor::AbstractMatrix{<:Union{RGB,RGBA}},
        modis_falsecolor::AbstractMatrix{<:Union{RGB,RGBA}},
        modis_cloud::AbstractMatrix{<:Union{RGB,RGBA}},
        labeled::AbstractMatrix,
        props::DataFrame,
        cloud_mask::AbstractMatrix,
        ice_mask::AbstractMatrix,
        landmask::AbstractMatrix,
        coastal_buffer_mask::AbstractMatrix,
        ift_version::VersionNumber = pkgversion(@__MODULE__),
        ift_archive_version::VersionNumber = VersionNumber("1.0.0"),
        reference::AbstractString = "https://doi.org/10.1016/j.rse.2019.111406",
        contact::AbstractString = "mmwilhelmus@brown.edu",
        creation_date::ZonedDateTime = now(tz"UTC"),
        ift_configuration::AbstractString = "",
    )

An object with results from a single segmentation to be saved as a proper
netCDF-4 file with [`save`](@ref).

Includes:

- References
  - `passtime`: the timepoint of the observation
  - `crs`: the CRS and geospatial data for the image, as returned by [`latlon`](@ref)
  - `modis_truecolor`: MODIS truecolor image (bands 1, 4, 3)
  - `modis_falsecolor`: MODIS falsecolor image (bands 7, 2, 1)
  - `modis_cloud`: the MODIS cloud RGB image
  - `ift_version`: the version of IceFloeTracker.jl used to save the file
  - `ift_archive_version`: the version of the file format (for this object, "1.0.0")
  - `reference`: a DOI for the dataset to which the file belongs
  - `contact`: contact information for the author
  - `creation_date`: ISO8601 timestamp of when the object was created (defaults to now)
  - `ift_configuration`: free-text description of the IFT configuration used to produce the file (e.g. the contents of `config.txt`)
- Images
  - `labeled`: the labeled image of connected components
  - `cloud_mask`: the cloud mask
  - `ice_mask`: the ice mask
  - `landmask`: the land mask
  - `coastal_buffer_mask`: the coastal buffer mask
- DataFrames
  - `props`: the measured properties of the floes

"""
@kwdef struct V1
    passtime::ZonedDateTime
    crs::NamedTuple
    modis_truecolor::AbstractMatrix{<:Union{RGB,RGBA}}
    modis_falsecolor::AbstractMatrix{<:Union{RGB,RGBA}}
    modis_cloud::AbstractMatrix{<:Union{RGB,RGBA}}
    labeled::AbstractMatrix
    props::DataFrame
    cloud_mask::AbstractMatrix
    ice_mask::AbstractMatrix
    landmask::AbstractMatrix
    coastal_buffer_mask::AbstractMatrix
    ift_version::VersionNumber = pkgversion(@__MODULE__)
    ift_archive_version::VersionNumber = VersionNumber("1.0.0")
    reference::AbstractString = "https://doi.org/10.1016/j.rse.2019.111406"
    contact::AbstractString = "mmwilhelmus@brown.edu"
    creation_date::ZonedDateTime = now(tz"UTC")
    ift_configuration::AbstractString = ""
end

"""
    save(path, V1(args...))

Write the [`V1`](@ref) object to storage as a proper netCDF-4 file.

The file uses `NCDatasets` and follows CF conventions. The `x`, `y`, and `time`
dimensions are defined at the root level and inherited by all groups.

Structure:

```
📦 netCDF-4 file
├─ 🏷️ ift_archive_version, ift_version, contact, reference   (global attributes)
├─ 🏷️ Description, label_variable, labeled_image     (floe-properties attributes)
├─ 🔢 geolocation   (scalar Int32, CRS grid-mapping variable)
│  └─ 🏷️ name, crs_wkt, spatial_ref, long_name, GeoTransform
├─ 🔢 x(x)          (projection x / easting coordinates, metres)
├─ 🔢 y(y)          (projection y / northing coordinates, metres)
├─ 🔢 time(time)    (seconds since 1970-01-01)
├─ 🔢 modis_falsecolor(x, y, band_modis_falsecolor)
│  └─ 🏷️ description, grid_mapping
├─ 🔢 modis_cloud(x, y, band_modis_cloud)
│  └─ 🏷️ description, grid_mapping
├─ 🔢 modis_truecolor(x, y, band_modis_truecolor)
│  └─ (same attributes as modis_falsecolor)
├─ 🔢 labeled_image(x, y)
│  └─ 🏷️ description, grid_mapping
├─ 🔢 cloud_mask(x, y)
├─ 🔢 coastal_buffer_mask(x, y)
├─ 🔢 ice_mask(x, y)
├─ 🔢 landmask(x, y)
├─ 🔢 floe_label(floe_label)   (links pixel values in labeled_image to rows here)
└─ 🔢 floe_<column>(floe_label) for every other column in the props DataFrame
```

"""
function save(output_path::AbstractString, s::V1;)
    ptsunix = Int64(Dates.datetime2unix(DateTime(s.passtime)))
    latlondata = s.crs

    crs_code = latlondata[:crs]
    crs_name = _get_crs_name(crs_code)
    projection_dataset_name = "geolocation"

    X = Float64.(latlondata[:X])
    Y = Float64.(latlondata[:Y])
    nx = length(X)
    ny = length(Y)

    NCDataset(output_path, "c") do ds
        # Global attributes
        ds.attrib["ift_archive_version"] = string(s.ift_archive_version)
        ds.attrib["ift_version"] = string(s.ift_version)
        ds.attrib["reference"] = s.reference
        ds.attrib["contact"] = s.contact
        ds.attrib["creation_date"] = Dates.format(
            s.creation_date, dateformat"yyyy-mm-ddTHH:MM:SSz"
        )
        ds.attrib["ift_configuration"] = s.ift_configuration

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
        vcrs.attrib["GeoTransform"] = join(latlondata[:geotransform], " ")
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
        nchannels_tc = size(channelview(s.modis_truecolor), 1)
        nchannels_fc = size(channelview(s.modis_falsecolor), 1)
        nchannels_mc = size(channelview(s.modis_cloud), 1)
        defDim(ds, "band_modis_truecolor", nchannels_tc)
        defDim(ds, "band_modis_falsecolor", nchannels_fc)
        defDim(ds, "band_modis_cloud", nchannels_mc)

        # band coordinate variables — string labels describing each channel
        _band_labels(names, n) =
            n == length(names) ? names : [names; fill("alpha", n - length(names))]
        tc_labels = _band_labels(TRUECOLOR_BANDLABELS, nchannels_tc)
        fc_labels = _band_labels(FALSECOLOR_BANDLABELS, nchannels_fc)
        vbtc = defVar(ds, "band_modis_truecolor", String, ("band_modis_truecolor",))
        vbtc.attrib["long_name"] = "MODIS truecolor band labels"
        vbtc[:] = tc_labels
        vbfc = defVar(ds, "band_modis_falsecolor", String, ("band_modis_falsecolor",))
        vbfc.attrib["long_name"] = "MODIS falsecolor band labels"
        vbfc[:] = fc_labels

        nc_create_color_dataset(
            ds,
            "modis_falsecolor",
            s.modis_falsecolor,
            "MODIS falsecolor image (bands 7, 2, 1)",
            projection_dataset_name,
            "band_modis_falsecolor",
        )
        nc_create_color_dataset(
            ds,
            "modis_truecolor",
            s.modis_truecolor,
            "MODIS truecolor image (bands 1, 4, 3)",
            projection_dataset_name,
            "band_modis_truecolor",
        )
        nc_create_color_dataset(
            ds,
            "modis_cloud",
            s.modis_cloud,
            "MODIS cloud RGB image",
            projection_dataset_name,
            "band_modis_cloud",
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
            "background cloud",
        )
        nc_create_mask_dataset(
            ds,
            "landmask",
            s.landmask,
            "Land mask. This mask is 1 for pixels classified as land, and 0 elsewhere.",
            projection_dataset_name,
            "background land",
        )
        nc_create_mask_dataset(
            ds,
            "coastal_buffer_mask",
            s.coastal_buffer_mask,
            "Coastal buffer mask. This mask is 1 for pixels within a specified " *
            "distance of the coast, and 0 elsewhere.",
            projection_dataset_name,
            "background coastal_buffer",
        )
        nc_create_mask_dataset(
            ds,
            "ice_mask",
            s.ice_mask,
            "Ice mask. This mask is 1 for pixels classified as ice, and 0 elsewhere.",
            projection_dataset_name,
            "background sea_ice",
        )

        # floe properties: one variable per property column
        nc_create_floe_properties(ds, s.props)
    end

    return nothing
end

"""
    nc_create_mask_dataset(grp, name, mask, description, projection_dataset_name, flag_meanings)

Define a `UInt8` binary mask variable with dimensions `(x, y)` in the NCDatasets
group `grp`. The parent group must supply the `x` and `y` dimensions.

CF `flag_values` (`[0x00, 0x01]`) and `flag_meanings` attributes are written to the
variable. `flag_meanings` should be a space-separated string of two labels, e.g.
`"background cloud"` or `"background sea_ice"`.
"""
function nc_create_mask_dataset(
    grp::NCDataset,
    name::AbstractString,
    mask::AbstractMatrix,
    description::AbstractString="",
    projection_dataset_name::AbstractString="geolocation",
    flag_meanings::AbstractString="background flag",
)
    mask_uint8 = permutedims(UInt8.(Bool.(mask)), (2, 1))  # (nx, ny)
    v = defVar(grp, name, UInt8, ("x", "y"))
    v.attrib["description"] = description
    v.attrib["grid_mapping"] = projection_dataset_name
    v.attrib["flag_values"] = UInt8[0, 1]
    v.attrib["flag_meanings"] = flag_meanings
    v[:, :] = mask_uint8
    return v
end

"""
    nc_create_labeled_dataset(grp, name, labeled, description, projection_dataset_name)

Define an integer-typed segmentation label variable with dimensions `(x, y)` in
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
    v.attrib["description"] = description
    v.attrib["grid_mapping"] = projection_dataset_name
    v[:, :, :] = img_raw
    return v
end

"""
    nc_create_floe_properties(grp, props)

Define one netCDF variable per column in the `props` DataFrame, all sharing a
`floe_label` dimension. Integer columns are stored as `Int64`; floating-point columns
as `Float64` with a `NaN` fill value; other columns as `String`. The `label`
variable links each floe back to pixel values in the labeled image via the group
attributes `label_variable` and `labeled_image`. CF `units`, `long_name`, and
`standard_name` attributes are attached to each recognized column.
"""
function nc_create_floe_properties(grp::NCDataset, props::DataFrame)
    props_ = convert_missing_to_nan(props)
    nfloes = nrow(props_)

    defDim(grp, "floe_label", nfloes)

    grp.attrib["label_variable"] = "floe_label"
    grp.attrib["labeled_image"] = "labeled_image"

    for col_name in names(props_)
        nc_name = "floe_" * col_name
        col_data = props_[!, col_name]
        T = eltype(col_data)

        v = _defVar(grp, nc_name, T, ("floe_label",))

        if haskey(COLUMN_ATTRIBUTES, col_name)
            for (k, val) in pairs(COLUMN_ATTRIBUTES[col_name])
                v.attrib[string(k)] = val
            end
        end

        if nfloes > 0
            v[:] = col_data
        end
    end
    return nothing
end

function _defVar(grp, name, ::Type{T}, dims) where {T<:Integer}
    defVar(grp, name, Int64, dims)
end
function _defVar(grp, name, ::Type{T}, dims) where {T<:AbstractFloat}
    defVar(grp, name, Float64, dims; fillvalue=NaN)
end
function _defVar(grp, name, ::Type, dims)
    defVar(grp, name, String, dims)
end

"""
    _load_v1(input_path)

Load a V1 netCDF-4 file written by [`save`](@ref) and return a [`V1`](@ref) object.
"""
function _load_v1(input_path::AbstractString)
    NCDataset(input_path, "r") do ds
        # Global attributes
        ift_version = VersionNumber(ds.attrib["ift_version"])
        ift_archive_version = VersionNumber(ds.attrib["ift_archive_version"])
        reference = ds.attrib["reference"]
        contact = ds.attrib["contact"]
        creation_date = ZonedDateTime(ds.attrib["creation_date"])
        ift_configuration = ds.attrib["ift_configuration"]

        # Reconstruct CRS data from the geolocation variable
        geoloc = ds["geolocation"]
        crs_wkt = geoloc.attrib["crs_wkt"]
        geotransform = parse.(Float64, split(geoloc.attrib["GeoTransform"]))
        crs = latlon(
            importWKT(crs_wkt),
            importEPSG(4326),
            geotransform,
            size(ds["y"].var)[1],
            size(ds["x"].var)[1],
        )

        # passtime — raw Int64 stored as "seconds since 1970-01-01"
        ptsunix = Int64(ds["time"].var[1])
        passtime = ZonedDateTime(Dates.unix2datetime(ptsunix), tz"UTC")

        # Colour images: stored as (nx, ny, nb) UInt8 → reconstruct Matrix{RGB/RGBA}
        function read_color(varname)
            raw = collect(permutedims(Array(ds[varname].var), (3, 2, 1)))  # (nb, ny, nx)
            nb = size(raw, 1)
            n0f8 = reinterpret(N0f8, raw)
            T = nb == 3 ? RGB{N0f8} : RGBA{N0f8}
            image = collect(colorview(T, n0f8))
            return image
        end
        modis_truecolor = read_color("modis_truecolor")
        modis_falsecolor = read_color("modis_falsecolor")
        modis_cloud = read_color("modis_cloud")

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

        # Floe properties: all variables sharing the "floe_label" dimension, in file order.
        # Strip the "floe_" prefix to recover original DataFrame column names.
        prop_cols = Pair{String,AbstractVector}[]
        for varname in keys(ds)
            v = ds[varname]
            "floe_label" ∈ dimnames(v) || continue
            col_name = startswith(varname, "floe_") ? varname[6:end] : varname
            # Use .var to bypass fill-value masking so NaN is preserved as NaN
            push!(prop_cols, col_name => Array(v.var))
        end
        props = isempty(prop_cols) ? DataFrame() : DataFrame(prop_cols)

        return V1(;
            passtime,
            crs,
            modis_truecolor,
            modis_falsecolor,
            modis_cloud,
            labeled,
            props,
            cloud_mask,
            ice_mask,
            landmask,
            coastal_buffer_mask,
            ift_version,
            ift_archive_version,
            reference,
            contact,
            creation_date,
            ift_configuration,
        )
    end
end

function _get_crs_name(crs_code)
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

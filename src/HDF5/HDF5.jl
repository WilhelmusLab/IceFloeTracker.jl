
module HDF5

export V1, load_hdf5, save_hdf5

using TimeZones, DataFrames

# These are implemented in ext/HDF5Ext.jl
function load_hdf5 end
function save_hdf5 end

"""
    IceFloeTracker.HDF5.V1(;
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
  - `truecolor_path`: the path to the truecolor image
  - `falsecolor_path`: the path to the falsecolor image
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
@kwdef struct V1
    passtime::ZonedDateTime
    crs_ref_image_path::AbstractString
    truecolor_path::AbstractString
    falsecolor_path::AbstractString
    labeled::AbstractMatrix
    props::DataFrame
    cloud_mask::AbstractMatrix
    ice_mask::AbstractMatrix
    landmask::AbstractMatrix
    coastal_buffer_mask::AbstractMatrix
    iftversion::VersionNumber = pkgversion(@__MODULE__)
    file_version::VersionNumber = VersionNumber("1.0.0")
    reference::AbstractString = "https://doi.org/10.1016/j.rse.2019.111406"
    contact::AbstractString = "mmwilhelmus@brown.edu"
end

end

import ArchGDAL:
    ISpatialRef,
    readraster,
    toPROJ4,
    importEPSG,
    toEPSG,
    getproj,
    getgeotransform,
    importWKT,
    toWKT
import Proj: Transformation

# px2xy function from creator of ArchGDAL https://github.com/yeesian/ArchGDAL.jl/issues/68
# Uses the geotransform from the image to map pixels to projection values
function px2xy(gt, x, y)
    xoff, a, b, yoff, d, e = gt
    xp = a * (x + 0.5) + b * (y + 0.5) + xoff
    yp = d * (x + 0.5) + e * (y + 0.5) + yoff
    return (xp, yp)
end

"""
    latlon(imgpath::AbstractString)

Reads the GeoTiff located at <imgpath>, extracts the coordinate reference system,
and produces a lookup table with for the column and row values in the same projection
as the GeoTiff and a 2D array for latitude and longitude.
"""
function latlon(imgpath::AbstractString)
    # rewrite of the latlon.py library with Julia
    # readraster automatically detects GeoTIFF
    im = readraster(imgpath)
    ref_ispatialref = importEPSG(4326)
    p_ispatialref = importWKT(getproj(im))
    nrows, ncols, _ = size(im)
    geotransform = getgeotransform(im)
    data = latlon(p_ispatialref, ref_ispatialref, geotransform, nrows, ncols)
    return data
end

function latlon(
    p_ispatialref::ISpatialRef,
    ref_ispatialref::ISpatialRef,
    geotransform::Union{NTuple{6,<:Real},Vector{<:Real}},
    nrows::Int,
    ncols::Int,
)
    ref = toPROJ4(ref_ispatialref)
    crs = toPROJ4(p_ispatialref)
    X = [px2xy(geotransform, i - 1, 0)[1] for i in 1:nrows]
    Y = [px2xy(geotransform, 0, j - 1)[2] for j in 1:ncols]
    trans = Transformation(crs, ref)
    lonlat = [trans(x, y) for y in Y, x in X]
    lon = first.(lonlat)
    lat = last.(lonlat)
    data = (;
        crs=toEPSG(p_ispatialref),
        crs_wkt=toWKT(p_ispatialref),
        longitude=lon,
        latitude=lat,
        X=X,
        Y=Y,
        geotransform=geotransform,
    )
    return data
end

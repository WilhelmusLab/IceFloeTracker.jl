import ArchGDAL as AG
using Proj

# px2xy function from creator of ArchGDAL https://github.com/yeesian/ArchGDAL.jl/issues/68
# Uses the geotransform from the image to map pixels to projection values
function px2xy(gt, x, y)
    xoff, a, b, yoff, d, e = gt
    xp = a * (x + 0.5) + b * (y + 0.5) + xoff
    yp = d * (x + 0.5) + e * (y + 0.5) + yoff
    (xp, yp)
end

"""
Reads the GeoTiff located at <imgpath>, extracts the coordinate reference system,
and produces a lookup table with for the column and row values in the same projection
as the GeoTiff and a 2D array for latitude and longitude.
"""
function latlon(imgpath)
# rewrite of the latlon.py library with Julia
    # readraster automatically detects GeoTIFF
    im = AG.readraster(imgpath) 

    # get a reference projection (WGS84)
    ref = AG.toPROJ4(AG.importEPSG(4326))

    # get the raster projection information and format it
    p = AG.getproj(im); 
    crs = AG.toPROJ4(AG.importWKT(p))

    nrows, ncols, _ = size(im)

    # Get the X and Y vectors
    gt = AG.getgeotransform(im)
    X = [px2xy(gt,i-1,0)[1] for i in 1:nrows]
    Y = [px2xy(gt,0,j-1)[2] for j in 1:ncols]

    # Similar to PyProj, going from source (polar stereo) to target (lat/lon)
    trans = Proj.Transformation(crs, ref)
    lonlat = [trans(x, y) for y in Y, x in X]
    lon = first.(lonlat)
    lat = last.(lonlat)

    return (crs=AG.toEPSG(AG.importWKT(p)),
            longitude=lon,
            latitude=lat,
            X=X,
            Y=Y)
end

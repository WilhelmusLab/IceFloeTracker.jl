import numpy as np
from pyproj import Transformer
import rasterio


def getlatlon(imgpath):
    """
    Get the longitude and latitude of the pixels in the image with path `imgpath`.
    :param imgpath: the path of the image
    """

    im = rasterio.open(imgpath)
    crs = im.crs.__str__()  # Coordinate reference system code
    nrows, ncols = im.shape
    rows, cols = np.meshgrid(np.arange(nrows), np.arange(ncols))
    xs, ys = rasterio.transform.xy(im.transform, rows, cols)
    # X and Y are the 1D index vectors
    X = np.array(xs)[:, 0]
    Y = np.array(ys)[0, :]
    ps_to_ll = Transformer.from_crs(im.crs, "WGS84", always_xy=True)
    lons, lats = ps_to_ll.transform(np.ravel(xs), np.ravel(ys))
    # longitude and latitude are 2D index arrays
    longitude = np.reshape(lons, (nrows, ncols))
    latitude = np.reshape(lats, (nrows, ncols))
    return {"crs": crs, "longitude": longitude, "latitude": latitude, "X": X, "Y": Y}

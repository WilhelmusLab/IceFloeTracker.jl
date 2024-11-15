import cv2 as cv
import numpy as np
from numpy.typing import ArrayLike


def imfill(binary_image: ArrayLike) -> ArrayLike:
    """Fill holes in a binary image.

    Args:
        binary_image: Input binary image.

    Returns:
        Output binary image with holes filled.
    """

    # Find contours with CV_RETR_CCOMP to handle holes (connected components)
    contours, hierarchy = cv.findContours(binary_image, cv.RETR_CCOMP, cv.CHAIN_APPROX_SIMPLE)

    # Create a blank image to draw the filled contours
    dst = np.zeros_like(binary_image)

    # Fill each contour, ensuring to handle holes correctly
    for i in range(len(contours)):
        if hierarchy[0][i][3] == -1:  # Draw only external contours
            cv.drawContours(dst, contours, i, (255, 255, 255), thickness=cv.FILLED)

    return dst
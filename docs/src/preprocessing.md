# Preprocessing

IFT operates on optical satellite imagery. The main functions are designed with "true color" and "false color" imagery in mind, and have thus far primarily been tested on imagery from the Moderate Resolution Imaging Spectroradiometer (MODIS) from the NASA _Aqua_ and _Terra_ satellites. The preprocessing routines mask land and cloud features, and aim to adjust and sharpen the remainder of the images to amplify the contrast along the edges of sea ice floes. We illustrate some of these tasks using the following real-world example, using MODIS images from the western Arctic Ocean. The functions use three different images, shown from left to right: a true color image, and a false color image, and a land mask. 

```@raw html
<img src="../assets/tc_fc_lm_example.png" width="600" alt="False color image, Lopez-Acosta cloud mask, Watkins cloud mask"/>
```

## Land masks

Landmask generation and dilation is handled by the function `create_landmask`. Landmask images from file are loaded as RGB matrices. This example uses an image from NASA EarthData landmask for Beaufort Sea.

```julia
using IceFloeTracker

rgb_landmask = IceFloeTracker.load("/path/to/landmask_image.tiff");
landmask_imgs = IceFloeTracker.create_landmask(rgb_landmask);
```

The `landmask_imgs` object includes a binary version of the original landmask and a dilated version, which helps to cover the complicated near-coastal regions.

```@raw html
<img src="../assets/landmask_example.png" width="400" alt="Landmask Example"/>
```

At the top, we have the original landmask TIFF, which has black and gray values. The middle image is the binary image, with land set to 0. At the bottom, we can see the dilated image using the default value of the structuring element. The default has radius 50, which results in a coastal mask of 12.5 km based on the 250 m pixel size of default MODIS imagery. Any structuring element compatible with Julia ImageMorphology is supported.

## Cloud masks
Clouds are near-ubiquitous in the summer Arctic Ocean. For segmenting sea ice imagery, it is generally required to identify cloud regions and process those separately from other regions. Currently IFT includes two closely related, customizable cloud mask functions which derive estimate cloud presence/absence from the false color images. 

Cloud masks are interpreted as binary images where `true` indicates the presence of cloud. They can be generated using the false color image and the `create_cloudmask()` function, with an cloud mask generating functor as an argument.

```julia
# Lopez-Acosta cloud mask
cloud_mask = create_cloudmask(falsecolor_image, LopezAcosta2019CloudMask())

# Watkins2025 cloud mask
cloud_mask = create_cloudmask(falsecolor_image, Watkins2025CloudMask())
```

The Lopez-Acosta 2019 cloud mask aims to only mask the brightest clouds, and is used in Pipelines where cloudy regions are further processed to enhance clarity. The Watkins 2025 cloud mask begins with the Lopez-Acosta 2019 algorithm, with different, stricter parameters, and includes morphological operations to remove speckle. This algorithm has lower tolerance for cloud cover.

```@raw html
<img src="../assets/cloudmask_example.png" width="400" alt="False color image, Lopez-Acosta cloud mask, Watkins cloud mask"/>
```
In the image, clouds are visible as bright patches in the false color image (top). The Lopez-Acosta 2019 cloud mask is in the middle, and the Watkins 2025 cloud mask is at the bottom. 

## Image filtering and adjustment
The IFT includes multiple functions for preparing an image for segmentation. Common tasks include equalizing the lighting in an image, smoothing noise away from boundaries, sharpening edges, and enhancing difference between floes and background ice. The primary image filters and image adjustment algorithms included in the `Filtering` module are

1. Nonlinear diffusion using the Perona-Malik algorithm. This algorithm performs diffusion using a heat equation weighted by a function image gradient. Thus, the diffusion is limited near edges and strong in object interiors. The algorithm was coded in Julia from the original Perona-Malik paper and includes both the inverse quadratic and exponential gradient functions.
```julia
pmd = PeronaMalikDiffusion(0.1, 0.1, 5, "exponential")
truecolor_diffused = nonlinear_diffusion(truecolor_image, pmd)
```
2. Adaptive histogram equalization. The Julia ImageContrastAdjustment library includes multiple methods for adjusting image histograms via the `adjust_histogram` function.
```julia
truecolor_equalized = adjust_histogram(truecolor_diffused,
                      AdaptiveEqualization(nbins=256, rblocks=8, cblocks=10, clip=0.8))
```
Note that the implementation of contrast limited adaptive histogram equalization differs from the version used in Matlab and in Python's `scikit-image` library, as the Julia version is based on a different source algorithm while having the same name. The Matlab and Python implementations are instead based on an algorithm published in (Graphics Gems IV)[https://github.com/erich666/GraphicsGems/blob/master/gemsiv/README]. We implemented the Graphic Gems version in Julia. The resulting function `ContrastLimitedAdaptiveHistogramEqualization` can be used as an alternative to `AdaptiveEqualization` in the `adjust_histogram` function.
3. Unsharp masking. This technique subtracts a Gaussian blurred copy of an image, resulting in sharper contrast at image object edges.
```julia
truecolor_sharpened =unsharp_mask(truecolor_equalized)
```

The image below illustrates the processing with a zoomed-in version of the example case. In the top row, we have the true color image on the left and the nonlinear diffusion image on the right. In the bottom row, we have the diffused and equalized image on the left and the diffused, equalized and sharpened image on the right.
```@raw html
<img src="../assets/filtering_example.png" width="600" alt="Original true color image, diffused image, equalized image, and sharpened image."/>
```
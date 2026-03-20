# Generate figures for the documentation website
using Pkg
HOME = "../."
Pkg.activate(HOME)
using IceFloeTracker
using Images

test_region = (1:2000, 1001:3000)

truecolor_image =
    float64.(
        RGB.(
            load(
                joinpath(
                    HOME,
                    "test/test_inputs/beaufort-chukchi-seas_truecolor.2020162.aqua.250m.tiff",
                ),
            )[test_region...]
        )
    )
falsecolor_image =
    float64.(
        RGB.(
            load(
                joinpath(
                    HOME,
                    "test/test_inputs/beaufort-chukchi-seas_falsecolor.2020162.aqua.250m.tiff",
                ),
            )[test_region...]
        )
    )

land_mask_img = load(joinpath(HOME, "test/test_inputs/landmask.tiff"))[test_region...]

save(
    joinpath(HOME, "docs/src/assets/tc_fc_lm_example.png"),
    imresize(
        mosaicview(
            truecolor_image,
            falsecolor_image,
            land_mask_img;
            nrow=1,
            npad=15,
            fillvalue=RGBA(0, 0, 0, 0),
        ),
        400,
        1206,
    ),
)

coastal_buffer, land_mask = create_landmask(land_mask_img)
# Cast to RGB or RGBA needed before fillvalue will work.
save(
    joinpath(HOME, "docs/src/assets/landmask_example.png"),
    imresize(
        mosaicview(
            land_mask_img,
            RGBA.(Gray.(land_mask)),
            RGBA.(Gray.(coastal_buffer));
            nrow=1,
            npad=15,
            fillvalue=RGBA(0, 0, 0, 0),
        ),
        400,
        1206,
    ),
)

cloudmask = create_cloudmask(falsecolor_image)
cm2 = create_cloudmask(falsecolor_image, Watkins2025CloudMask())
apply_landmask!(cm2, land_mask)
save(
    joinpath(HOME, "docs/src/assets/cloudmask_example.png"),
    imresize(
        mosaicview(
            falsecolor_image,
            RGBA.(Gray.(cloudmask)),
            RGBA.(Gray.(cm2));
            nrow=1,
            npad=15,
            fillvalue=RGBA(0, 0, 0, 0),
        ),
        400,
        1206,
    ),
)

# 1. Nonlinear diffusion using the Perona-Malik algorithm
pmd = PeronaMalikDiffusion(0.1, 0.1, 5, "exponential")
@time truecolor_diffused = nonlinear_diffusion(truecolor_image, pmd)

# 2. Adaptive histogram equalization
truecolor_equalized = adjust_histogram(
    truecolor_diffused,
    ContrastLimitedAdaptiveHistogramEqualization(;
        nbins=256, rblocks=10, cblocks=10, clip=10
    ),
)

# 3. Unsharp Mask and Convert to Grayscale
truecolor_sharpened = Gray.(unsharp_mask(truecolor_equalized))

# Apply masks
preprocessed_image = apply_landmask(truecolor_sharpened, land_mask)

zoom_region = (1001:2000, 501:1500)
save(
    joinpath(HOME, "docs/src/assets/filtering_example.png"),
    imresize(
        mosaicview(
            truecolor_diffused[zoom_region...],
            truecolor_equalized[zoom_region...],
            preprocessed_image[zoom_region...];
            nrow=1,
            npad=15,
            fillvalue=RGBA(0, 0, 0, 0),
        ),
        400,
        1206,
    ),
)

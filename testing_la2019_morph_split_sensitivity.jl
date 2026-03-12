using IceFloeTracker
using Images
using Random
using StatsBase
using DataFrames

function se_disk(r)
    se = [sum(abs.(c.I .- (r + 1)) .^ 2) for c in CartesianIndices((2*r + 1, 2*r + 1))]
    return sqrt.(se) .<= r
end

### Segmentation F Morph Split Function ####
# This function uses the hbreak/bridge, area opening, etc to split floes. Quite agressive.
# Here, I want to see how the different parameters change the output.
# The function exists in the FSPipeline module, but I reproduce it here in case I need to change details.
function morph_split_floes(binary_img, cloudmask; max_fill_area=1, min_area_opening=20, opening_strel=se_disk(4), bothat_strel=se_disk(1))
    leads_branched = hbreak(binary_img) |> branch
    leads_filled = .!imfill(.!leads_branched, 0:max_fill_area)
    leads_opened = branch(
        area_opening(leads_filled; min_area=min_area_opening, connectivity=2)
    )
    leads_bothat = bothat(leads_opened, bothat_strel) .> 0
    leads = convert(BitMatrix, (complement.(leads_bothat) .* leads_opened))
    area_opening!(leads, leads; min_area=min_area_opening, connectivity=2)

    floes = (fill_holes(leads) .* .!cloudmask) |> branch
    floes_opened = opening(floes, opening_strel)
    mreconstruct!(dilate, floes_opened, floes, floes_opened)
    return floes_opened
end

#### Load inputs ####
begin
    dataloc = "test/test_outputs/segmentation-IceFloeTracker.FSPipeline.Segment-006-baffin_bay-100km-20220530-terra-250m-2026-03-12-140258"
    validated = load(joinpath(dataloc, "validated_binary.png"))
    watershed_bdry = load(joinpath(dataloc, "watersheds_product.png"))
    seg_a = load(joinpath(dataloc, "segA.png"))
    ice_intersect = load(joinpath(dataloc, "segAB_intersect.png"))
    seg_f = load(joinpath(dataloc, "segF.png"))
    final = load(joinpath(dataloc, "final_floes.png"))
    preproc_gray = load(joinpath(dataloc, "sharpened_grayscale_image.png"))
    cloudmask = load(joinpath(dataloc, "cloud_mask.png"))
end

mosaicview(seg_a, ice_intersect, seg_f, watershed_bdry, nrow=1)

init_img = RGB.(seg_f .> 0) .+ RGB(0, 1, 0.3) .* Float64.(outlines)

# Illustrate steps with defaults
begin
    proc_images = []
    binary_img = deepcopy(ice_intersect)
    push!(proc_images, binary_img)
    max_fill_area=1
    min_area_opening=20
    opening_strel=se_disk(4)
    bothat_strel=strel_diamond((5,5))

    # 2
    leads_branched = hbreak(binary_img .> 0) |> branch
    push!(proc_images, Gray.(leads_branched))

    # 3
    leads_filled = .!imfill(.!leads_branched, 0:max_fill_area)
    push!(proc_images, Gray.(leads_filled))

    # 4
    leads_opened = branch(
        area_opening(leads_filled; min_area=min_area_opening, connectivity=2)
    )
    push!(proc_images, Gray.(leads_opened))

    # 5
    leads_bothat = bothat(leads_opened, bothat_strel) .> 0
    leads = convert(BitMatrix, (complement.(leads_bothat) .* leads_opened))
    push!(proc_images, Gray.(leads))

    # 6
    area_opening!(leads, leads; min_area=min_area_opening, connectivity=2)
    push!(proc_images, Gray.(leads))

    # 7
    floes = (fill_holes(leads) .* .!(cloudmask .> 0)) |> branch
    push!(proc_images, Gray.(floes))

    # 8
    floes_opened = opening(floes, opening_strel)
    push!(proc_images, Gray.(floes_opened))

    # 9
    mreconstruct!(dilate, floes_opened, floes, floes_opened)
    push!(proc_images, Gray.(floes_opened))
end

save("/Users/dwatkin2/Downloads/ift_split_stages.png", mosaicview([m for m in proc_images], nrow=3, rowmajor=true))

segF_overlay = RGB.(final) .+ RGB(1, 0, 0) .* outlines
segA_morphed_overlay = RGB.(floes_opened) .+ RGB(1, 0, 0) .* outlines;
mosaicview(segF_overlay, segA_morphed_overlay, nrow=1)



segA_overlay = 0.9 .* RGB.(seg_a) .+ RGB(1, 0, 0) .* outlines;
segA_overlay[outlines .> 0] .= RGB(1, 0, 0) 
segAB_overlay = 0.9 .* RGB.(ice_intersect) .+ RGB(1, 0, 0) .* outlines;
segAB_overlay[outlines .> 0] .= RGB(1, 0, 0) 
segF_overlay = 0.9 .* RGB.(seg_f) .+ RGB(1, 0, 0) .* outlines
segF_overlay[outlines .> 0] .= RGB(1, 0, 0) 
mosaicview(segA_overlay, segAB_overlay, segF_overlay, nrow=1)


segF_cleaned = opening(seg_f .> 0, se_disk(1))
segF_cleaned .= imfill(segF_cleaned, (0, 100)) # Clear small objects
segF_cleaned .= .!imfill(.!segF_cleaned, (0, 100)) # Fill small holes
mosaicview(seg_f, Gray.(segF_cleaned), nrow=1)

#### Sensitivity test: max fill area
begin
    outlines = IceFloeTracker.Morphology.bwperim(validated .> 0)
    images = [RGB.(float64.(morph_split_floes(seg_f .> 0, cloudmask .> 0; max_fill_area=idx))) for idx ∈ [1, 3, 5, 7, 9, 11]]
    for idx in range(1, length(images))
        im = images[idx]
        im = RGB.(float64.(im))
        im = im .+ (RGB(1, 0, 0) .* Float64.(outlines))
        images[idx] = im
    end
    mosaicview(images, nrow=2, rowmajor=true)
end

#### Sensitivity test: min area opening
begin
    outlines = IceFloeTracker.Morphology.bwperim(validated .> 0)
    images = [RGB.(float64.(morph_split_floes(seg_f .> 0, cloudmask .> 0; min_area_opening=idx)))
        for idx ∈ [5, 10, 15, 20, 25, 30]]
    for idx in range(1, length(images))
        im = images[idx]
        im = RGB.(float64.(im))
        im = im .+ (RGB(1, 0.5, 0) .* Float64.(outlines))
        images[idx] = im
    end
    mosaicview(images, nrow=2, rowmajor=true)
end
#### Sensitivity test: bottom hat transform
begin
    outlines = IceFloeTracker.Morphology.bwperim(validated .> 0)
    images = [RGB.(float64.(morph_split_floes(seg_f .> 0, cloudmask .> 0; bothat_strel=se_disk(idx))))
        for idx ∈ [2, 4, 6, 8, 10, 12]]
    for idx in range(1, length(images))
        im = images[idx]
        im = RGB.(float64.(im))
        im = im .+ (RGB(0.5, 0, 1) .* Float64.(outlines))
        images[idx] = im
    end
    mosaicview(images, nrow=2, rowmajor=true)
end









min_area_opening = 20
ice_leads = .!(watershed_bdry .> 0) .* (ice_intersect .> 0)
ice_leads .= .!area_opening(ice_leads; min_area=min_area_opening, connectivity=2)
ice_leads = opening(ice_leads, strel_diamond((3,3)))
# ice_leads[cloudmask .> 0] .= 1
function se_disk(r)
    se = [sum(abs.(c.I .- (r + 1)) .^ 2) for c in CartesianIndices((2*r + 1, 2*r + 1))]
    return sqrt.(se) .<= r
end

# Standard: dilated image as marker and intersect afterwards
begin
    marker_img = dilate(preproc_gray, strel_diamond((5, 5)))
    mreconstruct!(
        dilate, marker_img, complement.(marker_img), complement.(preproc_gray)
    )
    reconstructed_leads = complement.(marker_img) .* ice_leads
    overlay = RGB.(reconstructed_leads)
    outlines = validated_binary_floes(case) .- erode(validated_binary_floes(case));
    overlay[outlines .> 0] .= RGB(1,0,0)
    overlay
end

# Ice as markers?
begin
    marker_img = Gray.(complement.(ice_leads))
    mreconstruct!(
        dilate, marker_img, complement.(marker_img), complement.(preproc_gray)
    )
    reconstructed_leads = complement.(marker_img)

    overlay = RGB.(reconstructed_leads)
    outlines = validated_binary_floes(case) .- erode(validated_binary_floes(case));
    overlay[outlines .> 0] .= RGB(1,0,0)
    overlay
end

# Intersect with ice first?
begin
    masked_img = deepcopy(preproc_gray)
    masked_img[ice_leads .== 0] .= 0
    marker_img = dilate(masked_img, se_disk(1))
    mreconstruct!(
        dilate, marker_img, complement.(marker_img), complement.(masked_img)
    )
    reconstructed_leads = complement.(marker_img)

    overlay = RGB.(reconstructed_leads)
    outlines = validated_binary_floes(case) .- erode(validated_binary_floes(case));
    overlay[outlines .> 0] .= RGB(1,0,0)
    overlay
end

Gray.(ice_leads)
# alt - use the ice leads to change the image
marker_img1 = deepcopy(preproc_gray)
mreconstruct!(
    dilate, marker_img1, ice_leads, complement.(preproc_gray)
)
reconstructed_leads1 = complement.(marker_img1 .* ice_leads)

mosaicview(reconstructed_leads, reconstructed_leads1, nrow=1)
opened_img = opening(reconstructed_leads, FSPipeline.se_disk(3))
Gray.((Float64.(opened_img) .+ Float64.(reconstructed_leads)) ./ 2)

overlay = RGB.(opened_img)
overlay[outlines .> 0] .= RGB(1,0,0)
overlay


#######




@info "Building cloudmask"
# TODO: @hollandjg track down why the cloudmask is different for float32 vs float64 input images
    # dmw: I suspect it's likely it's roundoff error with the comparison to the ratio threshold, that's where I've seen differences
cloudmask = create_cloudmask(falsecolor_image)

# 2. Intermediate images
fc_masked = apply_landmask(falsecolor_image, coastal_buffer_mask)

segment = Watkins2026.Segment(; 
            preprocessing_function=Watkins2026.Preprocess(adapthisteq_params = (nbins=256, rblocks=2, cblocks=2, clip=0.99)),
            tile_size_pixels=400) # Could add flag for pixel size being too large


# 26 seconds on the first run, with a 400 pixel image
@time begin
    segments = segment(truecolor_image, falsecolor_image, coastal_buffer_mask)    
end

# view overlay
cview = view_seg_random(segments)
true_floes = validated_binary_floes(case) .> 0
outlines = true_floes .- erode(true_floes)
cview[segments.image_indexmap .== 0] .= RGB.(0,0,0)
cview[outlines .> 0] .= RGB.(1,1,1)

cview

la_segment = LopezAcosta2019Tiling.Segment()
@time begin
    la_segments = la_segment(truecolor_image, falsecolor_image, coastal_buffer_mask)    
end

##### Run on the comparison images in the Greenland Floe Scale Dataset folder #####

segment = Watkins2026.Segment(; 
            preprocessing_function=Watkins2026.Preprocess(adapthisteq_params = (nbins=256, rblocks=8, cblocks=4, clip=0.99)),
            tile_size_pixels=1000) # Could add flag for pixel size being too large

dataloc = "/Volumes/Research/ENG_Wilhelmus_Shared/group/IFT_fram_strait_dataset/"

case_info = DataFrame(
    year_folder = ["fram_strait-2003", "fram_strait-2006", "fram_strait-2009", "fram_strait-2013", "fram_strait-2014", "fram_strait-2018"],
    month_folder = ["fram_strait-20030331-20030501", "fram_strait-20060701-20060801", "fram_strait-20090331-20090501",
                    "fram_strait-20130331-20130501", "fram_strait-20140501-20140601", "fram_strait-20180601-20180701"],
    date = ["20030401", "20060710", "20090403", "20130424", "20140501", "20180607"],
    satellite = ["terra", "terra", "terra", "aqua", "aqua", "aqua"]
)
   
row = 1

# Load images
year_folder = case_info[row, "year_folder"]
month_folder = case_info[row, "month_folder"]
date = case_info[row, "date"]
satellite = case_info[row, "satellite"]

land_mask_img = Gray.(load(joinpath(dataloc, year_folder, month_folder, "landmask.tiff")))
true_color_img = RGB.(load(joinpath(dataloc, year_folder, month_folder, "truecolor",
                                 join([date, satellite, "truecolor", "250m", "tiff"], "."))))
false_color_img = RGB.(load(joinpath(dataloc, year_folder, month_folder, "falsecolor",
                                 join([date, satellite, "falsecolor", "250m", "tiff"], "."))))
comparison_img = Float64.(load(joinpath(dataloc, year_folder, month_folder, "labeled_raw",
                                 join([date, satellite, "labeled_raw", "250m", "tiff"], "."));));

land_mask = land_mask_img .> 0
coastal_buffer = dilate(land_mask, strel_box((51,51)));
comparison = SegmentedImage(true_color_img, label_components(comparison_img));

zoom = (1000:3000, 1000:3000)

@time begin
    new_segmentation = segment(true_color_img[zoom...], false_color_img[zoom...], land_mask[zoom...])    
end

# this took 730 seconds - yikes! 12 minutes for one image. Hope it did well!
saveloc = "/Users/dwatkin2/Documents/research/manuscripts/greenland_floe_scale_dataset/IFT_greenland_sea_dataset/figures/comparison/"*date*"/"
save(saveloc*date*"_watkins2026_tilesize1000_test.png", Gray.(new_segmentation.image_indexmap .> 0))
save(saveloc*date*"_watkins2026_tilesize1000_test_colorized.png", view_seg_random(new_segmentation))

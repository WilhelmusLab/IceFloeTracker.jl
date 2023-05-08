"""
Pair floes in the floe library using an equivalent implementation as in the MATLAB script `final_2020.m` from https://github.com/WilhelmusLab/ice-floe-tracker/blob/main/existing_code/final_2020.m.

# Arguments
- `indir`: path to directory containing the floe library (cropped floe masks for registration and correlation), extracted features from segmented images, and overpass elapsed times `delta_t`.
- `condition_thresholds`: 3-tuple of thresholds (each a named tuple) for deciding whether to match floe `i` from day `k` to floe j from day `k+1`.
- `mc_thresholds`: thresholds for area mismatch and psi-s shape correlation.

Following are the default set of thresholds `condition_thresholds` used for floe matching:
- Condition 1: time elapsed `dt` from image `k` to image `k+1` and distance between floes centroids `dist`: `t1=(dt = (30, 100, 1300), dist=(15, 30, 120))`

- Condition 2: area of floe i `area`, and the computed ratios for area, major axis,  minor axis, and convex area of floes `i` and `j` in days `k` and `k+1`, respectively: `t2=(area=1200, arearatio=0.28, majaxisratio=0.10, minaxisratio=0.12, convexarearatio=0.14)`

- Condition 3: as in the previous condition but set as follows: `
t3=(area=1200, arearatio=0.18, majaxisratio=0.07, minaxisratio=0.08, convexarearatio=0.09)`
"""
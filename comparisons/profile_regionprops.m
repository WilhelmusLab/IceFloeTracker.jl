% Generate areas for connected components for comparing with
% skimages.regionprops and profiling

test_image = imread("bw_new.tif");
stats = regionprops(test_image,'Area');
writematrix([stats.Area],'matlab_areas.csv')


ptCloud = pcread('body-v2.ply');
pcshow(ptCloud);

% STEP 1: construct graph

% STEP 2: sample graph (contour preserving sampling , random sampling , no subsampling)
sampling_density = 0.05;

% STEP 3: parallel segmentation (meanshift/quickshift - embarrassingly parallel):

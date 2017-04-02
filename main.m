% INITIALIZATION
addpath(genpath('libraries'),genpath('pcdata'));

%ptCloud = pcread('bun_zipper.ply');
% dataFile = fullfile(toolboxdir('vision'), 'visiondata', 'livingRoom.mat');
% load(dataFile);
% ptCloud = livingRoomData{44};

ptCloud = pcread('body-v2.ply');
gridStep = 20;
ptCloud = pcdownsample(ptCloud, 'gridAverage', gridStep);
pcshow(ptCloud);
%% PART 1 : GRAPH CONSTRUCTION:
%threshold = 1;
%disp('adding nodes to graph');
%G = graph_addnodes(ptCloud);
disp('adding edges to graph');
num_neighbors = 5;
sigma_sq = 10;
G = graph_addedges(ptCloud,num_neighbors,sigma_sq);
Graph_pc = graph_addnodeattributes(G,ptCloud);
%% PART 2 : GRAPH RESAMPLING:


% STEP 2: sample graph (contour preserving sampling , random sampling , no subsampling)
sampling_density = 0.05;










%% PART 3 : SEGMENTATION:


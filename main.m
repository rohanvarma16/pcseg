% INITIALIZATION
addpath(genpath('libraries'),genpath('pcdata'));

%ptCloud = pcread('bun_zipper.ply');
dataFile = fullfile(toolboxdir('vision'), 'visiondata', 'livingRoom.mat');
load(dataFile);
ptCloud = livingRoomData{1};


% flatten and remove NaN's:

%ptCloud_i = clean_pc(ptCloud);
ptCloud_i = load('ptCloud_room.mat');
ptCloud_i = ptCloud_i.ptCloud;
iscolor = 1;
% ptCloud = pcread('body-v2.ply');
% gridStep = 10;
% iscolor = 1;
% ptCloud = pcdownsample(ptCloud, 'gridAverage', gridStep);
%ptCloud = pcread('teapot.ply');
figure(1);
pcshow(ptCloud_i);
gridStep = 0.02;
ptCloud = pcdownsample(ptCloud_i, 'gridAverage', gridStep);
figure(2)
pcshow(ptCloud);
%% PART 1 : GRAPH CONSTRUCTION:
%threshold = 1;
%disp('adding nodes to graph');
%G = graph_addnodes(ptCloud);
num_neighbors = 10;
sigma_sq = 0.01;
disp('adding edges to graph'); 
G = graph_addedges(ptCloud,num_neighbors,sigma_sq);
disp('adding node attributes graph'); 
Graph_pc = graph_addnodeattributes(G,ptCloud,iscolor);
%% PART 2 : GRAPH RESAMPLING:
Graph_pc_imp = graph_calcimp(Graph_pc, ptCloud);
sampling_density = 0.5;
sample_index = randsample(ptCloud.Count, round(sampling_density * ptCloud.Count), true, Graph_pc_imp.Nodes.Imp);
%construct new point cloud from sample_index:
xyz_rs = ptCloud.Location(sample_index,:);
if(iscolor)
    color_rs = ptCloud.Color(sample_index,:);
    ptCloud_rs = pointCloud(xyz_rs, 'Color',color_rs);
else
    ptCloud_rs = pointCloud(xyz_rs);
end
figure(1);
pcshow(ptCloud);
figure(2);
pcshow(ptCloud_rs);



%% PART 3 : SEGMENTATION









%% PART 3 : SEGMENTATION:


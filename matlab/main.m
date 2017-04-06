% INITIALIZATION
 addpath(genpath('libraries'),genpath('pcdata'));
% 
% %ptCloud = pcread('bun_zipper.ply');
dataFile = fullfile(toolboxdir('vision'), 'visiondata', 'livingRoom.mat');
load(dataFile);
ptCloud = livingRoomData{1};
 iscolor = 1;


%flatten and remove NaN's:

ptCloud_i = clean_pc(ptCloud);

ptCloud_i = load('ptCloud_room.mat');
ptCloud_i = ptCloud_i.ptCloud;
% 
% ptCloud_i = pcread('body-v2.ply');
% figure(1); 
%  pcshow(ptCloud_i);
% %for body
%  gridStep = 10;
% iscolor = 1;
% ptCloud = pcdownsample(ptCloud_i, 'gridAverage', gridStep);

%for room:
  gridStep = 0.02;
 ptCloud = pcdownsample(ptCloud_i, 'gridAverage', gridStep);
iscolor = 1;
figure(2)
pcshow(ptCloud_i);
%% PART 1 : GRAPH CONSTRUCTION:

num_neighbors = 20;
%for room
sigma_sq = 0.01;

%for body
%sigma_sq = 10;


disp('adding edges to graph'); 
G = graph_addedges(ptCloud,num_neighbors,sigma_sq);
disp('adding node attributes graph'); 
Graph_pc = graph_addnodeattributes(G,ptCloud,iscolor);
%% PART 2 : GRAPH RESAMPLING:
Graph_pc_imp = graph_calcimp(Graph_pc, ptCloud);
sampling_density = 0.2;
sample_index = randsample(numnodes(Graph_pc_imp), round(sampling_density * numnodes(Graph_pc_imp)), true, Graph_pc_imp.Nodes.Imp);
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

G_rs = subgraph(Graph_pc_imp, unique(sample_index));
%% plot top K points

% np = 60000;
% [sort_wt, sort_idx] = sort(Graph_pc_imp.Nodes.Imp,'descend');
% sample_index = sort_idx(1:np);
% xyz_rs = ptCloud.Location(sample_index,:);
% if(iscolor)
%     color_rs = ptCloud.Color(sample_index,:);
%     ptCloud_rs = pointCloud(xyz_rs, 'Color',color_rs);
% else
%     ptCloud_rs = pointCloud(xyz_rs);
% end
% figure(1);
% pcshow(ptCloud);
% figure(2);
% pcshow(ptCloud_rs);


%% PART 3 : SEGMENTATION
[P_n, parent_x, d_x] = segmentation(G_rs);







%% PART 3 : SEGMENTATION:


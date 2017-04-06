%% point cloud segmentation given a point cloud/graph
function [P_n, parent_x, d_x] = segmentation(G)
numpts = numnodes(G);
P_n = zeros(numpts,1);

%% compute density:
for i = 1:numpts
    disp(i);
    nodeid = G.Nodes.Name(i);
    nodeid = num2str(nodeid{1});
    neighbors_i = neighbors(G, nodeid);
    sum = 0;
    for j = 1:length(neighbors_i)
        weight = G.Edges.Weight(findedge(G,nodeid, neighbors_i{j}));
        sum = sum + weight;
    end
    P_n(i) = sum;
end

%% link neighbors:
%tau = 0;
d_x = ones(numpts,1)*inf;
parent_x = zeros(numpts,1);

for i = 1:numpts
    disp(i);
    nodeid = G.Nodes.Name(i);
    nodeid = num2str(nodeid{1});   
    neighbors_i = neighbors(G, nodeid);
    min_distance = inf;
    nodenames = G.Nodes.Name;
    for j = 1:length(neighbors_i)
        nbr = neighbors_i{j};
        nbr_ind = find( strcmp(nodenames,num2str(nbr)) == 1);
        weight = G.Edges.Dist(findedge(G,nodeid, neighbors_i{j}));
       if (P_n(nbr_ind) > P_n(i) && weight < min_distance)
           d_x(i) = weight;
           parent_x(i) = nbr_ind;
           min_distance = weight;
       end
    end
end

end

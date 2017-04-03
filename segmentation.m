%% point cloud segmentation given a point cloud/graph
P_n = zeros(ptCloud.Count,1);

%% compute density:
for i = 1:ptCloud.Count
    disp(i);
    neighbors_i = neighbors(G, int2str(i));
    sum = 0;
    for j = 1:length(neighbors_i)
        weight = G.Edges.Weight(findedge(G,int2str(i), neighbors_i{j}));
        sum = sum + weight;
    end
    P_n(i) = sum;
end

%% link neighbors:
tau = 0;
d_x = zeros(ptCloud.Count,1);
parent_x = zeros(ptCloud.Count, 1);

for i = 1:ptCloud.Count
    disp(i);
    neighbors_i = neighbors(G, int2str(i));
    min_distance = inf;
    for j = 1:length(neighbors_i)
        weight = G.Edges.Weight(findedge(G,int2str(i), neighbors_i{j}));
       if (P_n(j) > P_n(i) && weight < min_distance)
           d_x(i) = weight;
           parent_x(i) = j;
           min_distance = weight;
       end
    end
end

function G = graph_addedges(ptCloud, num_neighbors, sigma_sq)
G = graph();
for i = 1:ptCloud.Count
   G = addnode(G,int2str(i)); 
end

for i = 1:ptCloud.Count
    disp(i);
    [neighborhood, distances] = findNearestNeighbors(ptCloud, ...
        [ptCloud.Location(i,1),ptCloud.Location(i,2), ptCloud.Location(i,3)], ...
         num_neighbors, 'MaxLeafChecks',100);
     T = table();
    for j = 1:length(distances)
        if(distances(j) > 0)
            if(findedge(G,int2str(i),int2str(neighborhood(j))) == 0)
                weight = exp(-1*distances(j)/sigma_sq);
                disp(weight);
                edgeprops = table( [ {int2str(i)}, {int2str(neighborhood(j))}], ...
                    weight, 'VariableNames', {'EndNodes', 'Weight'} );
                T = [T ; edgeprops];
            end
        end
    end
    if(size(T,1) > 0)
        G = addedge(G,T);
    end
end
end
